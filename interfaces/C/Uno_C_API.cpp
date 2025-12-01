// Copyright (c) 2025 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cstring>
#include <iostream>
#include <streambuf>
#include "Uno_C_API.h"
#include "../UserModel.hpp"
#include "Uno.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "model/Model.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Presets.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/CollectionAdapter.hpp"
#include "symbolic/Range.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"
#include "tools/NumberModelEvaluations.hpp"
#include "tools/UserCallbacks.hpp"

using namespace uno;

using CUserModel = UserModel<Objective, ObjectiveGradient, Constraints, Jacobian, JacobianOperator, JacobianTransposedOperator,
   Hessian, HessianOperator, std::vector<double>, void*>;

// UnoModel contains an instance of UserModel and complies with the Model interface
class UnoModel: public Model {
public:
   explicit UnoModel(const CUserModel& user_model):
         Model("C model", static_cast<size_t>(user_model.number_variables), static_cast<size_t>(user_model.number_constraints),
            static_cast<double>(user_model.optimization_sense), user_model.lagrangian_sign_convention),
         user_model(user_model),
         nonlinear_constraints(this->number_constraints),
         equality_constraints_collection(this->equality_constraints),
         inequality_constraints_collection(this->inequality_constraints) {
      this->find_fixed_variables(this->fixed_variables);
      this->partition_constraints(this->equality_constraints, this->inequality_constraints);
   }

   [[nodiscard]] ProblemType get_problem_type() const override {
      return this->user_model.problem_type;
   }

   // availability of linear operators
   [[nodiscard]] bool has_jacobian_operator() const override {
      return (this->user_model.jacobian_operator != nullptr);
   }

   [[nodiscard]] bool has_jacobian_transposed_operator() const override {
      return (this->user_model.jacobian_transposed_operator != nullptr);
   }

   [[nodiscard]] bool has_hessian_operator() const override {
      return (this->user_model.lagrangian_hessian_operator != nullptr);
   }

   [[nodiscard]] bool has_hessian_matrix() const override {
      return (this->user_model.lagrangian_hessian != nullptr);
   }

   // function evaluations
   [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override {
      double objective_value{0.};
      if (this->user_model.objective_function != nullptr) {
         const uno_int return_code = this->user_model.objective_function(this->user_model.number_variables, x.data(),
            &objective_value, this->user_model.user_data);
         if (0 < return_code) {
            throw FunctionEvaluationError();
         }
         objective_value *= this->optimization_sense;
         ++this->number_model_evaluations.objective;
      }
      return objective_value;
   }

   void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override {
      if (this->user_model.constraint_functions != nullptr) {
         const uno_int return_code = this->user_model.constraint_functions(this->user_model.number_variables,
            this->user_model.number_constraints, x.data(), constraints.data(), this->user_model.user_data);
         if (0 < return_code) {
            throw FunctionEvaluationError();
         }
         ++this->number_model_evaluations.constraints;
      }
   }

   // dense objective gradient
   void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
      if (this->user_model.objective_gradient != nullptr) {
         const uno_int return_code = this->user_model.objective_gradient(this->user_model.number_variables, x.data(),
            gradient.data(), this->user_model.user_data);
         if (0 < return_code) {
            throw GradientEvaluationError();
         }
         for (size_t variable_index: Range(this->number_variables)) {
            gradient[variable_index] *= this->optimization_sense;
         }
         ++this->number_model_evaluations.objective_gradient;
      }
   }

   // sparsity patterns of Jacobian and Hessian
   void compute_constraint_jacobian_sparsity(uno_int * row_indices, uno_int * column_indices, uno_int solver_indexing,
                                             MatrixOrder /*matrix_order*/) const override {
      // copy the indices of the user sparsity patterns to the Uno vectors
      for (size_t nonzero_index: Range(static_cast<size_t>(this->user_model.number_jacobian_nonzeros))) {
         row_indices[nonzero_index] = this->user_model.jacobian_row_indices[nonzero_index];
         column_indices[nonzero_index] = this->user_model.jacobian_column_indices[nonzero_index];
      }
      // TODO matrix_order

      // handle the solver indexing
      if (this->user_model.base_indexing != solver_indexing) {
         const uno_int indexing_difference = solver_indexing - this->user_model.base_indexing;
         for (size_t nonzero_index: Range(static_cast<size_t>(this->user_model.number_jacobian_nonzeros))) {
            row_indices[nonzero_index] += indexing_difference;
            column_indices[nonzero_index] += indexing_difference;
         }
      }
   }

   void compute_hessian_sparsity(uno_int *row_indices, uno_int *column_indices, uno_int solver_indexing) const override {
      // copy the indices of the user sparsity patterns to the Uno vectors
      const size_t number_hessian_nonzeros = this->number_hessian_nonzeros();
      for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
         row_indices[nonzero_index] = this->user_model.hessian_row_indices[nonzero_index];
         column_indices[nonzero_index] = this->user_model.hessian_column_indices[nonzero_index];
      }

      // handle the solver indexing
      if (this->user_model.base_indexing != solver_indexing) {
         const uno_int indexing_difference = solver_indexing - this->user_model.base_indexing;
         for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
            row_indices[nonzero_index] += indexing_difference;
            column_indices[nonzero_index] += indexing_difference;
         }
      }
   }

   // numerical evaluations of Jacobian and Hessian
   void evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const override {
      if (this->user_model.constraint_jacobian != nullptr) {
         const uno_int return_code = this->user_model.constraint_jacobian(this->user_model.number_variables,
            this->user_model.number_jacobian_nonzeros, x.data(), jacobian_values, this->user_model.user_data);
         if (0 < return_code) {
            throw GradientEvaluationError();
         }
         ++this->number_model_evaluations.jacobian;
      }
   }

   void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         double* hessian_values) const override {
      if (this->user_model.lagrangian_hessian != nullptr) {
         objective_multiplier *= this->optimization_sense;
         // if the model has a different sign convention for the Lagrangian than Uno, flip the signs of the multipliers
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            const_cast<Vector<double>&>(multipliers).scale(-1.);
         }
         const uno_int number_hessian_nonzeros = static_cast<uno_int>(this->number_hessian_nonzeros());
         const uno_int return_code = this->user_model.lagrangian_hessian(this->user_model.number_variables,
            this->user_model.number_constraints, number_hessian_nonzeros, x.data(), objective_multiplier,
            multipliers.data(), hessian_values, this->user_model.user_data);
         // flip the signs of the multipliers back
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            const_cast<Vector<double>&>(multipliers).scale(-1.);
         }
         if (0 < return_code) {
            throw HessianEvaluationError();
         }
         ++this->number_model_evaluations.hessian;
      }
      else {
         throw std::runtime_error("evaluate_lagrangian_hessian not implemented");
      }
   }

   void compute_jacobian_vector_product(const double* x, const double* vector, double* result) const override {
      if (this->user_model.jacobian_operator != nullptr) {
         const uno_int return_code = this->user_model.jacobian_operator(this->user_model.number_variables,
            this->user_model.number_constraints, x, true, vector, result, this->user_model.user_data);
         if (0 < return_code) {
            throw GradientEvaluationError();
         }
      }
      else {
         throw std::runtime_error("compute_jacobian_vector_product not implemented");
      }
   }

   void compute_jacobian_transposed_vector_product(const double* x, const double* vector, double* result) const override {
      if (this->user_model.jacobian_transposed_operator != nullptr) {
         const uno_int return_code = this->user_model.jacobian_transposed_operator(this->user_model.number_variables,
            this->user_model.number_constraints, x, true, vector, result, this->user_model.user_data);
         if (0 < return_code) {
            throw GradientEvaluationError();
         }
      }
      else {
         throw std::runtime_error("compute_jacobian_transposed_vector_product not implemented");
      }
   }

   void compute_hessian_vector_product(const double* x, const double* vector, double objective_multiplier, const Vector<double>& multipliers,
         double* result) const override {
      if (this->user_model.lagrangian_hessian_operator != nullptr) {
         objective_multiplier *= this->optimization_sense;
         // if the model has a different sign convention for the Lagrangian than Uno, flip the signs of the multipliers
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            const_cast<Vector<double>&>(multipliers).scale(-1.);
         }
         const uno_int return_code = this->user_model.lagrangian_hessian_operator(this->user_model.number_variables,
            this->user_model.number_constraints, x, true, objective_multiplier, multipliers.data(), vector, result, this->user_model.user_data);
         // flip the signs of the multipliers back
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            const_cast<Vector<double>&>(multipliers).scale(-1.);
         }
         if (0 < return_code) {
            throw HessianEvaluationError();
         }
      }
      else {
         throw std::runtime_error("compute_hessian_vector_product not implemented");
      }
   }

   [[nodiscard]] double variable_lower_bound(size_t variable_index) const override {
      return this->user_model.variables_lower_bounds[variable_index];
   }

   [[nodiscard]] double variable_upper_bound(size_t variable_index) const override {
      return this->user_model.variables_upper_bounds[variable_index];
   }

   [[nodiscard]] const SparseVector<size_t>& get_slacks() const override {
      return this->slacks;
   }

   [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override {
      return this->fixed_variables;
   }

   [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override {
      return this->user_model.constraints_lower_bounds[constraint_index];
   }

   [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override {
      return this->user_model.constraints_upper_bounds[constraint_index];
   }

   [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override {
      return this->equality_constraints_collection;
   }

   [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override {
      return this->inequality_constraints_collection;
   }

   [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override {
      return this->linear_constraints;
   }

   [[nodiscard]] const Collection<size_t>& get_nonlinear_constraints() const override {
      return this->nonlinear_constraints;
   }

   void initial_primal_point(Vector<double>& x) const override {
      // copy the initial primal point
      for (size_t variable_index: Range(static_cast<size_t>(this->user_model.number_variables))) {
         x[variable_index] = this->user_model.initial_primal_iterate[variable_index];
      }
   }

   void initial_dual_point(Vector<double>& multipliers) const override {
      // copy the initial dual point
      for (size_t constraint_index: Range(static_cast<size_t>(this->user_model.number_constraints))) {
         multipliers[constraint_index] = this->user_model.initial_dual_iterate[constraint_index];
      }
      if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
         multipliers.scale(-1.);
      }
   }

   void postprocess_solution(Iterate& /*iterate*/) const override {
      // do nothing
   }

   [[nodiscard]] size_t number_jacobian_nonzeros() const override {
      return static_cast<size_t>(this->user_model.number_jacobian_nonzeros);
   }

   [[nodiscard]] size_t number_hessian_nonzeros() const override {
      if (this->user_model.number_hessian_nonzeros.has_value()) {
         return static_cast<size_t>(*this->user_model.number_hessian_nonzeros);
      }
      else {
         throw std::runtime_error("The number of Hessian nonzeros is not available in UnoModel");
      }
   }

   [[nodiscard]] size_t number_model_objective_evaluations() const override {
      return this->number_model_evaluations.objective;
   }

   [[nodiscard]] size_t number_model_constraints_evaluations() const override {
      return this->number_model_evaluations.constraints;
   }

   [[nodiscard]] size_t number_model_objective_gradient_evaluations() const override {
      return this->number_model_evaluations.objective_gradient;
   }

   [[nodiscard]] size_t number_model_jacobian_evaluations() const override {
      return this->number_model_evaluations.jacobian;
   }

   [[nodiscard]] size_t number_model_hessian_evaluations() const override {
      return this->number_model_evaluations.hessian;
   }

   void reset_number_evaluations() const override {
      this->number_model_evaluations.reset();
   }

protected:
   const CUserModel& user_model;
   mutable NumberModelEvaluations number_model_evaluations{};
   const SparseVector<size_t> slacks{};
   Vector<size_t> fixed_variables{};
   const ForwardRange linear_constraints{0};
   const ForwardRange nonlinear_constraints;
   std::vector<size_t> equality_constraints;
   CollectionAdapter<std::vector<size_t>> equality_constraints_collection;
   std::vector<size_t> inequality_constraints;
   CollectionAdapter<std::vector<size_t>> inequality_constraints_collection;
};

class CUserCallbacks: public UserCallbacks {
public:
   CUserCallbacks(NotifyAcceptableIterateUserCallback notify_acceptable_iterate_callback,
      TerminationUserCallback user_termination_callback, void* user_data): UserCallbacks(),
         notify_acceptable_iterate_callback(notify_acceptable_iterate_callback),
         user_termination_callback(user_termination_callback),
         user_data(user_data) { };

   void notify_acceptable_iterate(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier,
         double primal_feasibility_residual, double stationarity_residual, double complementarity_residual) override {
      if (this->notify_acceptable_iterate_callback != nullptr) {
         this->notify_acceptable_iterate_callback(static_cast<uno_int>(primals.size()),
            static_cast<uno_int>(multipliers.constraints.size()), primals.data(), multipliers.lower_bounds.data(),
            multipliers.upper_bounds.data(), multipliers.constraints.data(), objective_multiplier, primal_feasibility_residual,
            stationarity_residual, complementarity_residual, this->user_data);
      }
   }

   bool user_termination(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier,
         double primal_feasibility_residual, double stationarity_residual, double complementarity_residual) override {
      if (this->user_termination_callback) {
         return this->user_termination_callback(static_cast<uno_int>(primals.size()),
            static_cast<uno_int>(multipliers.constraints.size()), primals.data(), multipliers.lower_bounds.data(),
            multipliers.upper_bounds.data(), multipliers.constraints.data(), objective_multiplier, primal_feasibility_residual,
            stationarity_residual, complementarity_residual, this->user_data);
      }
      else {
         return false; // never terminate
      }
   }

private:
   NotifyAcceptableIterateUserCallback notify_acceptable_iterate_callback;
   TerminationUserCallback user_termination_callback;
   void* user_data;
};

// std::streambuf wrapper around LoggerStreamUserCallback
class CStreamBuffer : public std::streambuf {
public:
   explicit CStreamBuffer(LoggerStreamUserCallback logger_stream_callback, void* user_data, std::size_t buffer_size) :
         logger_stream_callback(logger_stream_callback), user_data(user_data) {
      // allocate output buffer and set stream buffer pointer
      this->buffer = new char[buffer_size];
      this->setp(this->buffer, this->buffer + buffer_size - 1);
   }
   ~CStreamBuffer() override {
      // flush remaining data and release buffer memory
      this->sync();
      delete[] this->buffer;
   }

protected:
   // called on buffer overflow
   int overflow(int character) override {
      if (character != EOF) {
            // insert the character into the buffer
            *this->pptr() = traits_type::to_char_type(character);
            this->pbump(1);
      }
      // return EOF for error
      return (this->flush_buffer() == 0) ? character : EOF;
   }

   int sync() override {
      return this->flush_buffer();
   }

private:
   LoggerStreamUserCallback logger_stream_callback;
   void* user_data;
   char* buffer;

   // flush buffer to the logger callback
   int flush_buffer() {
      // check for invalid stream callback
      if (!this->logger_stream_callback) {
         return -1;
      }
      std::ptrdiff_t current_used_buffer_size = this->pptr() - this->pbase();
      if (current_used_buffer_size > 0) {
         // call user logger callback
         const uno_int callback_result = this->logger_stream_callback(this->pbase(), static_cast<uno_int>(current_used_buffer_size),
            this->user_data);
         if (callback_result != static_cast<uno_int>(current_used_buffer_size)) {
            return -1;
         }
         // move buffer pointer
         this->pbump(static_cast<int>(-current_used_buffer_size));
      }
      return 0;
   }
};

// std::ostream wrapper around LoggerStreamUserCallback
class COStream : public std::ostream {
public:
   COStream(LoggerStreamUserCallback logger_stream_callback, void* user_data, std::size_t buffer_size = 1024) : // 1024 default buffer size
      std::ostream(&this->buffer), buffer(logger_stream_callback, user_data, buffer_size) { }

private:
   // internal stream buffer that sends output to the LoggerStreamUserCallback
   CStreamBuffer buffer;
};

COStream* c_ostream = nullptr;

struct Solver {
   Uno* solver;
   Options* options;
   UserCallbacks* user_callbacks;
   Result* result;
};

void uno_get_version(uno_int* major, uno_int* minor, uno_int* patch) {
   *major = UNO_VERSION_MAJOR;
   *minor = UNO_VERSION_MINOR;
   *patch = UNO_VERSION_PATCH;
}

void* uno_create_model(const char* problem_type, uno_int number_variables, const double* variables_lower_bounds,
      const double* variables_upper_bounds, uno_int base_indexing) {
   if (number_variables <= 0) {
      WARNING << "Please specify a positive number of variables."  << std::endl;
      return nullptr;
   }
   if (strcmp(problem_type, UNO_PROBLEM_LINEAR) != 0 && strcmp(problem_type, UNO_PROBLEM_QUADRATIC) != 0 &&
         strcmp(problem_type, UNO_PROBLEM_NONLINEAR) != 0) {
      WARNING << "Please specify a valid problem type."  << std::endl;
      return nullptr;
   }
   if (base_indexing != UNO_ZERO_BASED_INDEXING && base_indexing != UNO_ONE_BASED_INDEXING) {
      WARNING << "Please specify a valid base indexing." << std::endl;
      return nullptr;
   }
   CUserModel* user_model = new CUserModel(problem_type, number_variables, base_indexing);
   // copy the bounds internally
   const size_t unsigned_number_variables = static_cast<size_t>(number_variables);
   user_model->variables_lower_bounds.resize(unsigned_number_variables);
   user_model->variables_upper_bounds.resize(unsigned_number_variables);
   for (size_t variable_index: Range(unsigned_number_variables)) {
      user_model->variables_lower_bounds[variable_index] = (variables_lower_bounds != nullptr) ?
         variables_lower_bounds[variable_index] : -INF<double>;
      user_model->variables_upper_bounds[variable_index] = (variables_upper_bounds != nullptr) ?
         variables_upper_bounds[variable_index] : INF<double>;
   }
   // create the initial primal point
   user_model->initial_primal_iterate.resize(unsigned_number_variables);
   for (size_t variable_index: Range(unsigned_number_variables)) {
      user_model->initial_primal_iterate[variable_index] = 0.;
   }
   return user_model;
}

bool uno_set_objective(void* model, uno_int optimization_sense, Objective objective_function,
      ObjectiveGradient objective_gradient) {
   if (optimization_sense != UNO_MINIMIZE && optimization_sense != UNO_MAXIMIZE) {
      WARNING << "Please specify a valid objective sense."  << std::endl;
      return false;
   }
   if (model == nullptr) {
      WARNING << "Please specify a valid model."  << std::endl;
      return false;
   }
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->optimization_sense = optimization_sense;
   user_model->objective_function = objective_function;
   user_model->objective_gradient = objective_gradient;
   return true;
}

bool uno_set_constraints(void* model, uno_int number_constraints, Constraints constraint_functions,
      const double* constraints_lower_bounds, const double* constraints_upper_bounds, uno_int number_jacobian_nonzeros,
      const uno_int* jacobian_row_indices, const uno_int* jacobian_column_indices, Jacobian constraint_jacobian) {
   if (number_constraints <= 0) {
      WARNING << "Please specify a positive number of constraints."  << std::endl;
      return false;
   }
   if (model == nullptr) {
      WARNING << "Please specify a valid model."  << std::endl;
      return false;
   }
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->number_constraints = number_constraints;
   user_model->constraint_functions = constraint_functions;
   // copy the bounds internally
   const size_t unsigned_number_constraints = static_cast<size_t>(number_constraints);
   user_model->constraints_lower_bounds.resize(unsigned_number_constraints);
   user_model->constraints_upper_bounds.resize(unsigned_number_constraints);
   for (size_t constraint_index: Range(unsigned_number_constraints)) {
      user_model->constraints_lower_bounds[constraint_index] = (constraints_lower_bounds != nullptr) ?
         constraints_lower_bounds[constraint_index] : -INF<double>;
      user_model->constraints_upper_bounds[constraint_index] = (constraints_upper_bounds != nullptr) ?
         constraints_upper_bounds[constraint_index] : INF<double>;
   }
   user_model->number_jacobian_nonzeros = number_jacobian_nonzeros;
   // copy the Jacobian sparsity to allow the calling code to dispose of its vectors
   user_model->jacobian_row_indices.resize(static_cast<size_t>(number_jacobian_nonzeros));
   user_model->jacobian_column_indices.resize(static_cast<size_t>(number_jacobian_nonzeros));
   for (size_t index: Range(static_cast<size_t>(number_jacobian_nonzeros))) {
      user_model->jacobian_row_indices[index] = jacobian_row_indices[index];
      user_model->jacobian_column_indices[index] = jacobian_column_indices[index];
   }
   user_model->constraint_jacobian = constraint_jacobian;
   // create the initial dual point
   user_model->initial_dual_iterate.resize(unsigned_number_constraints);
   for (size_t constraint_index: Range(unsigned_number_constraints)) {
      user_model->initial_dual_iterate[constraint_index] = 0.;
   }
   return true;
}

bool uno_set_jacobian_operator(void* model, JacobianOperator jacobian_operator) {
   if (model == nullptr) {
      WARNING << "Please specify a valid model."  << std::endl;
      return false;
   }
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->jacobian_operator = jacobian_operator;
   return true;
}

bool uno_set_jacobian_transposed_operator(void* model, JacobianTransposedOperator jacobian_transposed_operator) {
   if (model == nullptr) {
      WARNING << "Please specify a valid model."  << std::endl;
      return false;
   }
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->jacobian_transposed_operator = jacobian_transposed_operator;
   return true;
}

bool uno_set_lagrangian_hessian(void* model, uno_int number_hessian_nonzeros, char hessian_triangular_part,
      const uno_int* hessian_row_indices, const uno_int* hessian_column_indices, Hessian lagrangian_hessian,
      double lagrangian_sign_convention) {
   if (number_hessian_nonzeros <= 0) {
      WARNING << "Please specify a positive number of Lagrangian Hessian nonzeros."  << std::endl;
      return false;
   }
   if (hessian_triangular_part != UNO_LOWER_TRIANGLE && hessian_triangular_part != UNO_UPPER_TRIANGLE) {
      WARNING << "Please specify a correct Hessian triangle in {'" << UNO_LOWER_TRIANGLE << "', '" <<
         UNO_UPPER_TRIANGLE << "'}."  << std::endl;
      return false;
   }
   if (lagrangian_sign_convention != UNO_MULTIPLIER_NEGATIVE && lagrangian_sign_convention != UNO_MULTIPLIER_POSITIVE) {
      WARNING << "Please specify a correct Lagrangian sign convention in {" << UNO_MULTIPLIER_NEGATIVE << ", " <<
         UNO_MULTIPLIER_POSITIVE << "}."  << std::endl;
      return false;
   }
   if (model == nullptr) {
      WARNING << "Please specify a valid model."  << std::endl;
      return false;
   }
   CUserModel* user_model = static_cast<CUserModel*>(model);
   // make sure that the sign convention is consistent with that of the Hessian operator
   if (user_model->lagrangian_hessian_operator != nullptr && user_model->lagrangian_sign_convention != lagrangian_sign_convention) {
      WARNING << "Please specify a Lagrangian sign convention consistent with that of the Hessian operator."  << std::endl;
      return false;
   }
   user_model->number_hessian_nonzeros = number_hessian_nonzeros;
   // copy the Hessian sparsity to allow the calling code to dispose of its vectors
   // from now on, we only maintain the lower triangle
   user_model->hessian_row_indices.resize(static_cast<size_t>(number_hessian_nonzeros));
   user_model->hessian_column_indices.resize(static_cast<size_t>(number_hessian_nonzeros));
   const bool lower_triangle = (hessian_triangular_part == UNO_LOWER_TRIANGLE);
   for (size_t index: Range(static_cast<size_t>(number_hessian_nonzeros))) {
      user_model->hessian_row_indices[index] = lower_triangle ? hessian_row_indices[index] : hessian_column_indices[index];
      user_model->hessian_column_indices[index] = lower_triangle ? hessian_column_indices[index] : hessian_row_indices[index];
   }
   user_model->hessian_triangular_part = UNO_LOWER_TRIANGLE;
   user_model->lagrangian_hessian = lagrangian_hessian;
   user_model->lagrangian_sign_convention = lagrangian_sign_convention;
   return true;
}

bool uno_set_lagrangian_hessian_operator(void* model, HessianOperator lagrangian_hessian_operator, double lagrangian_sign_convention) {
   if (lagrangian_sign_convention != UNO_MULTIPLIER_NEGATIVE && lagrangian_sign_convention != UNO_MULTIPLIER_POSITIVE) {
      WARNING << "Please specify a Lagrangian sign convention in {" << UNO_MULTIPLIER_NEGATIVE << ", " <<
         UNO_MULTIPLIER_POSITIVE << "}."  << std::endl;
      return false;
   }
   if (model == nullptr) {
      WARNING << "Please specify a valid model."  << std::endl;
      return false;
   }
   CUserModel* user_model = static_cast<CUserModel*>(model);
   // make sure that the sign convention is consistent with that of the Hessian function
   if (user_model->lagrangian_hessian != nullptr && user_model->lagrangian_sign_convention != lagrangian_sign_convention) {
      WARNING << "Please specify a Lagrangian sign convention consistent with that of the Hessian function."  << std::endl;
      return false;
   }
   user_model->lagrangian_hessian_operator = lagrangian_hessian_operator;
   user_model->lagrangian_sign_convention = lagrangian_sign_convention;
   return true;
}

bool uno_set_user_data(void* model, void* user_data) {
   if (model == nullptr) {
      WARNING << "Please specify a valid model."  << std::endl;
      return false;
   }
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->user_data = user_data;
   return true;
}

bool uno_set_initial_primal_iterate_component(void* model, uno_int index, double initial_primal_component) {
   if (model == nullptr) {
      WARNING << "Please specify a valid model."  << std::endl;
      return false;
   }
   CUserModel* user_model = static_cast<CUserModel*>(model);
   if (0 <= index && index < user_model->number_variables) {
      const size_t unsigned_index = static_cast<size_t>(index);
      user_model->initial_primal_iterate[unsigned_index] = initial_primal_component;
      return true;
   }
   return false;
}

bool uno_set_initial_dual_iterate_component(void* model, uno_int index, double initial_dual_component) {
   if (model == nullptr) {
      WARNING << "Please specify a valid model."  << std::endl;
      return false;
   }
   CUserModel* user_model = static_cast<CUserModel*>(model);
   if (0 <= index && index < user_model->number_constraints) {
      const size_t unsigned_index = static_cast<size_t>(index);
      user_model->initial_dual_iterate[unsigned_index] = initial_dual_component;
      return true;
   }
   return false;
}

bool uno_set_initial_primal_iterate(void* model, const double* initial_primal_iterate) {
   if (model == nullptr) {
      WARNING << "Please specify a valid model."  << std::endl;
      return false;
   }
   if (initial_primal_iterate != nullptr) {
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // copy the initial primal point
      for (size_t variable_index: Range(static_cast<size_t>(user_model->number_variables))) {
         user_model->initial_primal_iterate[variable_index] = initial_primal_iterate[variable_index];
      }
      return true;
   }
   return false;
}

bool uno_set_initial_dual_iterate(void* model, const double* initial_dual_iterate) {
   if (model == nullptr) {
      WARNING << "Please specify a valid model."  << std::endl;
      return false;
   }
   if (initial_dual_iterate != nullptr) {
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // copy the initial dual point
      for (size_t constraint_index: Range(static_cast<size_t>(user_model->number_constraints))) {
         user_model->initial_dual_iterate[constraint_index] = initial_dual_iterate[constraint_index];
      }
      return true;
   }
   return false;
}

void* uno_create_solver() {
   // default options
   Options* options = new Options;
   DefaultOptions::load(*options);
   // add default preset
   const Options preset_options = Presets::get_preset_options(std::nullopt);
   options->overwrite_with(preset_options);

   // default user callbacks
   UserCallbacks* user_callbacks = new NoUserCallbacks;

   // Uno solver
   Uno* uno_solver = new Uno;
   Solver* solver = new Solver{uno_solver, options, user_callbacks, nullptr}; // no result yet
   return solver;
}

bool uno_set_solver_integer_option(void* solver, const char* option_name, uno_int option_value) {
   if (solver == nullptr) {
      WARNING << "Please specify a valid solver."  << std::endl;
      return false;
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   uno_solver->options->set_integer(option_name, option_value);
   return true;
}

bool uno_set_solver_double_option(void* solver, const char* option_name, double option_value) {
   if (solver == nullptr) {
      WARNING << "Please specify a valid solver."  << std::endl;
      return false;
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   uno_solver->options->set_double(option_name, option_value);
   return true;
}

bool uno_set_solver_bool_option(void* solver, const char* option_name, bool option_value) {
   if (solver == nullptr) {
      WARNING << "Please specify a valid solver."  << std::endl;
      return false;
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   uno_solver->options->set_bool(option_name, option_value);
   return true;
}

bool uno_set_solver_string_option(void* solver, const char* option_name, const char* option_value) {
   if (solver == nullptr) {
      WARNING << "Please specify a valid solver."  << std::endl;
      return false;
   }
   // handle the preset separately
   if (strcmp(option_name, "preset") == 0) {
      uno_set_solver_preset(solver, option_value);
   }
   else {
      Solver* uno_solver = static_cast<Solver*>(solver);
      uno_solver->options->set_string(option_name, option_value);
   }
   return true;
}

uno_int uno_get_solver_option_type(void* solver, const char* option_name) {
   if (solver == nullptr) {
      WARNING << "Please specify a valid solver."  << std::endl;
      return false;
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   try {
      return static_cast<uno_int>(uno_solver->options->get_option_type(option_name));
   }
   catch(const std::out_of_range&) {
      return UNO_OPTION_TYPE_NOT_FOUND;
   }
}

bool uno_load_solver_option_file(void* solver, const char* file_name) {
   if (solver == nullptr) {
      WARNING << "Please specify a valid solver."  << std::endl;
      return false;
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   uno_solver->options->overwrite_with(uno::Options::load_option_file(file_name));
   return true;
}

bool uno_set_solver_preset(void* solver, const char* preset_name) {
   if (solver == nullptr) {
      WARNING << "Please specify a valid solver."  << std::endl;
      return false;
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   Presets::set(*uno_solver->options, preset_name);
   return true;
}

bool uno_set_solver_callbacks(void* solver, NotifyAcceptableIterateUserCallback notify_acceptable_iterate_callback,
      TerminationUserCallback user_termination_callback, void* user_data) {
   if (solver == nullptr) {
      WARNING << "Please specify a valid solver."  << std::endl;
      return false;
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   delete uno_solver->user_callbacks; // delete the previous callbacks
   uno_solver->user_callbacks = new CUserCallbacks(notify_acceptable_iterate_callback, user_termination_callback, user_data);
   return true;
}

bool uno_set_logger_stream_callback(LoggerStreamUserCallback logger_stream_callback, void* user_data) {
   delete c_ostream;
   c_ostream = new COStream(logger_stream_callback, user_data);
   Logger::set_stream(*c_ostream);
   return true;
}

bool uno_reset_logger_stream() {
   delete c_ostream;
   c_ostream = nullptr;
   Logger::set_stream(std::cout);
   return true;
}

void uno_optimize(void* solver, void* model) {
   if (model == nullptr) {
      throw std::runtime_error("Please specify a valid model.");
   }
   if (solver == nullptr) {
      throw std::runtime_error("Please specify a valid solver.");
   }
   CUserModel* user_model = static_cast<CUserModel*>(model);
   if (!user_model->objective_function && !user_model->constraint_functions) {
      WARNING << "Please specify at least an objective or constraints."  << std::endl;
      return;
   }
   Solver* uno_solver = static_cast<Solver*>(solver);

   // create an instance of UnoModel, a subclass of Model, and solve the model using Uno
   const UnoModel uno_model(*user_model);
   Logger::set_logger(uno_solver->options->get_string("logger"));
   Result result = uno_solver->solver->solve(uno_model, *uno_solver->options, *uno_solver->user_callbacks);
   // clean up the previous result (if any)
   delete uno_solver->result;
   // move the new result into uno_solver
   uno_solver->result = new Result(std::move(result));
   // flush the logger
   Logger::flush();
}

double uno_get_solver_double_option(void* solver, const char* option_name) {
   if (solver == nullptr) {
      throw std::runtime_error("Please specify a valid solver.");
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   return uno_solver->options->get_double(option_name);
}

int uno_get_solver_integer_option(void* solver, const char* option_name) {
   if (solver == nullptr) {
      throw std::runtime_error("Please specify a valid solver.");
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   return uno_solver->options->get_int(option_name);
}

bool uno_get_solver_bool_option(void* solver, const char* option_name) {
   if (solver == nullptr) {
      throw std::runtime_error("Please specify a valid solver.");
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   return uno_solver->options->get_bool(option_name);
}

const char* uno_get_solver_string_option(void* solver, const char* option_name) {
   if (solver == nullptr) {
      throw std::runtime_error("Please specify a valid solver.");
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   return uno_solver->options->get_string(option_name).c_str();
}

// auxiliary function
Result* uno_get_result(void* solver) {
   if (solver == nullptr) {
      throw std::runtime_error("Please specify a valid solver.");
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   if (uno_solver->result == nullptr) {
      throw std::runtime_error("The result is not available.");
   }
   return uno_solver->result;
}

uno_int uno_get_optimization_status(void* solver) {
   const Result* result = uno_get_result(solver);
   return static_cast<uno_int>(result->optimization_status);
}

uno_int uno_get_solution_status(void* solver) {
   const Result* result = uno_get_result(solver);
   return static_cast<uno_int>(result->solution_status);
}

double uno_get_solution_objective(void* solver) {
   const Result* result = uno_get_result(solver);
   return result->solution_objective;
}

double uno_get_primal_solution_component(void* solver, uno_int index) {
   const Result* result = uno_get_result(solver);
   const size_t unsigned_index = static_cast<size_t>(index);
   if (index < 0 || result->number_variables <= unsigned_index) {
      throw std::runtime_error("The index is not valid.");
   }
   return result->primal_solution[unsigned_index];
}

double uno_get_constraint_dual_solution_component(void* solver, uno_int index) {
   const Result* result = uno_get_result(solver);
   const size_t unsigned_index = static_cast<size_t>(index);
   if (index < 0 || result->number_constraints <= unsigned_index) {
      throw std::runtime_error("The index is not valid.");
   }
   return result->constraint_dual_solution[unsigned_index];
}

double uno_get_lower_bound_dual_solution_component(void* solver, uno_int index) {
   const Result* result = uno_get_result(solver);
   const size_t unsigned_index = static_cast<size_t>(index);
   if (index < 0 || result->number_variables <= unsigned_index) {
      throw std::runtime_error("The index is not valid.");
   }
   return result->lower_bound_dual_solution[unsigned_index];
}

double uno_get_upper_bound_dual_solution_component(void* solver, uno_int index) {
   const Result* result = uno_get_result(solver);
   const size_t unsigned_index = static_cast<size_t>(index);
   if (index < 0 || result->number_variables <= unsigned_index) {
      throw std::runtime_error("The index is not valid.");
   }
   return result->upper_bound_dual_solution[unsigned_index];
}

void uno_get_primal_solution(void* solver, double* primal_solution) {
   const Result* result = uno_get_result(solver);
   for (size_t variable_index: Range(result->number_variables)) {
      primal_solution[variable_index] = result->primal_solution[variable_index];
   }
}

void uno_get_constraint_dual_solution(void* solver, double* constraint_dual_solution) {
   const Result* result = uno_get_result(solver);
   for (size_t constraint_index: Range(result->number_constraints)) {
      constraint_dual_solution[constraint_index] = result->constraint_dual_solution[constraint_index];
   }
}

void uno_get_lower_bound_dual_solution(void* solver, double* lower_bound_dual_solution) {
   const Result* result = uno_get_result(solver);
   for (size_t variable_index: Range(result->number_variables)) {
      lower_bound_dual_solution[variable_index] = result->lower_bound_dual_solution[variable_index];
   }
}

void uno_get_upper_bound_dual_solution(void* solver, double* upper_bound_dual_solution) {
   const Result* result = uno_get_result(solver);
   for (size_t variable_index: Range(result->number_variables)) {
      upper_bound_dual_solution[variable_index] = result->upper_bound_dual_solution[variable_index];
   }
}

double uno_get_solution_primal_feasibility(void* solver) {
   const Result* result = uno_get_result(solver);
   return result->solution_primal_feasibility;
}

double uno_get_solution_stationarity(void* solver) {
   const Result* result = uno_get_result(solver);
   return result->solution_stationarity;
}

double uno_get_solution_complementarity(void* solver) {
   const Result* result = uno_get_result(solver);
   return result->solution_complementarity;
}

uno_int uno_get_number_iterations(void* solver) {
   const Result* result = uno_get_result(solver);
   return static_cast<uno_int>(result->number_iterations);
}

double uno_get_cpu_time(void* solver) {
   const Result* result = uno_get_result(solver);
   return result->cpu_time;
}

uno_int uno_get_number_objective_evaluations(void* solver) {
   const Result* result = uno_get_result(solver);
   return static_cast<uno_int>(result->number_objective_evaluations);
}

uno_int uno_get_number_constraint_evaluations(void* solver) {
   const Result* result = uno_get_result(solver);
   return static_cast<uno_int>(result->number_constraint_evaluations);
}

uno_int uno_get_number_objective_gradient_evaluations(void* solver) {
   const Result* result = uno_get_result(solver);
   return static_cast<uno_int>(result->number_objective_gradient_evaluations);
}

uno_int uno_get_number_jacobian_evaluations(void* solver) {
   const Result* result = uno_get_result(solver);
   return static_cast<uno_int>(result->number_jacobian_evaluations);
}

uno_int uno_get_number_hessian_evaluations(void* solver) {
   const Result* result = uno_get_result(solver);
   return static_cast<uno_int>(result->number_hessian_evaluations);
}

uno_int uno_get_number_subproblem_solved_evaluations(void* solver) {
   const Result* result = uno_get_result(solver);
   return static_cast<uno_int>(result->number_subproblems_solved);
}

void uno_destroy_model(void* model) {
   if (model != nullptr) {
      delete static_cast<CUserModel*>(model);
   }
}

void uno_destroy_solver(void* solver) {
   if (solver != nullptr) {
      Solver* uno_solver = static_cast<Solver*>(solver);
      delete uno_solver->solver;
      delete uno_solver->options;
      if (uno_solver->result != nullptr) {
         delete uno_solver->result;
      }
      delete uno_solver;
   }
}