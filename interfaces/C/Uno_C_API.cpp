// Copyright (c) 2025 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <memory>
#include <new>
#include <streambuf>
#include "Uno_C_API.h"
#include "../UserModel.hpp"
#include "Uno.hpp"
#include "linear_algebra/Vector.hpp"
#include "model/Model.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/IntegerRange.hpp"
#include "symbolic/Range.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"
#include "tools/NumberModelEvaluations.hpp"
#include "tools/UserCallbacks.hpp"

using namespace uno;

using OptionIterator = decltype(Options::option_types)::const_iterator;

using CUserModel = UserModel<uno_objective_callback, uno_objective_gradient_callback, uno_constraints_callback,
   uno_constraints_jacobian_callback, uno_constraints_jacobian_operator_callback, uno_constraints_jacobian_transposed_operator_callback,
   uno_lagrangian_hessian_callback, uno_lagrangian_hessian_operator_callback, std::vector<double>, void*>;

// UnoModel contains an instance of UserModel and complies with the Model interface
class UnoModel: public Model {
public:
   explicit UnoModel(const CUserModel& user_model):
         Model(user_model.name.empty() ? "C model" : user_model.name,
            static_cast<size_t>(user_model.number_variables), static_cast<size_t>(user_model.number_constraints),
            static_cast<double>(user_model.optimization_sense), static_cast<double>(user_model.lagrangian_sign_convention),
            user_model.base_indexing),
         user_model(user_model),
         nonlinear_constraints(this->number_constraints),
         variables_lower_bounds(normalize_infinite_bounds(user_model.variables_lower_bounds)),
         variables_upper_bounds(normalize_infinite_bounds(user_model.variables_upper_bounds)),
         constraints_lower_bounds(normalize_infinite_bounds(user_model.constraints_lower_bounds)),
         constraints_upper_bounds(normalize_infinite_bounds(user_model.constraints_upper_bounds)) {
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
   View<const uno_int> get_jacobian_row_indices() const override {
      return view(this->user_model.jacobian_row_indices.data(), this->number_jacobian_nonzeros());
   }

   View<const uno_int> get_jacobian_column_indices() const override {
      return view(this->user_model.jacobian_column_indices.data(), this->number_jacobian_nonzeros());
   }

   void compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override {
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
   void evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const override {
      if (this->user_model.jacobian != nullptr) {
         const uno_int return_code = this->user_model.jacobian(this->user_model.number_variables,
            this->user_model.number_jacobian_nonzeros, x.data(), jacobian_values, this->user_model.user_data);
         if (0 < return_code) {
            throw JacobianEvaluationError();
         }
         ++this->number_model_evaluations.jacobian;
      }
   }

   void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         View<double> hessian_values) const override {
      if (this->user_model.lagrangian_hessian != nullptr) {
         objective_multiplier *= this->optimization_sense;
         // if the model has a different sign convention for the Lagrangian than Uno, flip the signs of the multipliers
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            const_cast<Vector<double>&>(multipliers).scale(-1.);
         }
         const uno_int number_hessian_nonzeros = static_cast<uno_int>(this->number_hessian_nonzeros());
         const uno_int return_code = this->user_model.lagrangian_hessian(this->user_model.number_variables,
            this->user_model.number_constraints, number_hessian_nonzeros, x.data(), objective_multiplier,
            multipliers.data(), hessian_values.data(), this->user_model.user_data);
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
            throw JacobianEvaluationError();
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
            throw JacobianEvaluationError();
         }
      }
      else {
         throw std::runtime_error("compute_jacobian_transposed_vector_product not implemented");
      }
   }

   void compute_hessian_vector_product(View<const double> x, View<const double> vector, double objective_multiplier,
         const Vector<double>& multipliers, View<double> result) const override {
      if (this->user_model.lagrangian_hessian_operator != nullptr) {
         objective_multiplier *= this->optimization_sense;
         // if the model has a different sign convention for the Lagrangian than Uno, flip the signs of the multipliers
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            const_cast<Vector<double>&>(multipliers).scale(-1.);
         }
         const uno_int return_code = this->user_model.lagrangian_hessian_operator(this->user_model.number_variables,
            this->user_model.number_constraints, x.data(), true, objective_multiplier, multipliers.data(), vector.data(),
            result.data(), this->user_model.user_data);
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

   [[nodiscard]] const std::vector<double>& get_variables_lower_bounds() const override {
      return this->variables_lower_bounds;
   }

   [[nodiscard]] const std::vector<double>& get_variables_upper_bounds() const override {
      return this->variables_upper_bounds;
   }

   [[nodiscard]] const std::vector<double>& get_constraints_lower_bounds() const override {
      return this->constraints_lower_bounds;
   }

   [[nodiscard]] const std::vector<double>& get_constraints_upper_bounds() const override {
      return this->constraints_upper_bounds;
   }

   [[nodiscard]] const std::vector<size_t>& get_equality_constraints() const override {
      return this->equality_constraints;
   }

   [[nodiscard]] const std::vector<size_t>& get_inequality_constraints() const override {
      return this->inequality_constraints;
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

   void postprocess_solution(Iterate& /*iterate*/, Evaluations& /*evaluations*/) const override {
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
   const IntegerRange linear_constraints{0};
   const IntegerRange nonlinear_constraints;
   std::vector<size_t> equality_constraints;
   std::vector<size_t> inequality_constraints;
   // user bounds with infinite values made explicit
   const std::vector<double> variables_lower_bounds;
   const std::vector<double> variables_upper_bounds;
   const std::vector<double> constraints_lower_bounds;
   const std::vector<double> constraints_upper_bounds;
};

class CUserCallbacks: public UserCallbacks {
public:
   CUserCallbacks(uno_notify_acceptable_iterate_callback notify_acceptable_iterate_callback,
      uno_termination_callback termination_callback, void* user_data): UserCallbacks(),
         notify_acceptable_iterate_callback(notify_acceptable_iterate_callback),
         termination_callback(termination_callback),
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

   bool termination(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier,
         double primal_feasibility_residual, double stationarity_residual, double complementarity_residual) override {
      if (this->termination_callback) {
         return this->termination_callback(static_cast<uno_int>(primals.size()),
            static_cast<uno_int>(multipliers.constraints.size()), primals.data(), multipliers.lower_bounds.data(),
            multipliers.upper_bounds.data(), multipliers.constraints.data(), objective_multiplier, primal_feasibility_residual,
            stationarity_residual, complementarity_residual, this->user_data) != 0;
      }
      else {
         return false; // no user termination
      }
   }

private:
   uno_notify_acceptable_iterate_callback notify_acceptable_iterate_callback;
   uno_termination_callback termination_callback;
   void* user_data;
};

// std::streambuf wrapper around uno_logger_stream_callback
class CStreamBuffer : public std::streambuf {
public:
   explicit CStreamBuffer(uno_logger_stream_callback logger_stream_callback, void* user_data, std::size_t buffer_size) :
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
   uno_logger_stream_callback logger_stream_callback;
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

// std::ostream wrapper around uno_logger_stream_callback
class COStream : public std::ostream {
public:
   COStream(uno_logger_stream_callback logger_stream_callback, void* user_data, std::size_t buffer_size = 1024) : // 1024 default buffer size
      std::ostream(&this->buffer), buffer(logger_stream_callback, user_data, buffer_size) { }

private:
   // internal stream buffer that sends output to the uno_logger_stream_callback
   CStreamBuffer buffer;
};

COStream* c_ostream = nullptr;

struct Solver {
   std::unique_ptr<Uno> solver{std::make_unique<Uno>()};
   Options user_options{};
   Options preset_options{}; // empty for "auto"
   Options default_options{};
   std::unique_ptr<UserCallbacks> user_callbacks{std::make_unique<NoUserCallbacks>()};
   std::unique_ptr<Result> result{}; // the result of the last optimization, if it succeeded
   bool options_consistent_with_solve{false}; // is true only after the solve, false whenever new options are set
   // the options used by the last solve (defaults -> preset, including the values picked by the "auto" preset -> user).
   // Valid only if options_consistent_with_solve
   Options solve_options{};
};

namespace {
   // reports an error without throwing
   void report_error(const char* message) noexcept {
      try {
         WARNING << message << std::endl;
      }
      catch (...) {
      }
   }

   // runs a function and returns its value, or "failure_value" if it throws: no exception crosses the C interface
   template <typename Value, typename Function>
   Value guard(Value failure_value, Function&& function) noexcept {
      try {
         return function();
      }
      catch (const std::exception& exception) {
         report_error(exception.what());
      }
      catch (...) {
         report_error("An unknown error occurred");
      }
      return failure_value;
   }
} // namespace

void uno_get_version(uno_int* major, uno_int* minor, uno_int* patch) UNO_NOEXCEPT {
   if (major != nullptr) {
      *major = UNO_VERSION_MAJOR;
   }
   if (minor != nullptr) {
      *minor = UNO_VERSION_MINOR;
   }
   if (patch != nullptr) {
      *patch = UNO_VERSION_PATCH;
   }
}

void* uno_create_model(const char* problem_type, uno_int number_variables, const double* variables_lower_bounds,
      const double* variables_upper_bounds, uno_int base_indexing) UNO_NOEXCEPT {
   return guard(static_cast<void*>(nullptr), [&]() -> void* {
      void* model = uno_create_unconstrained_model(problem_type, number_variables, base_indexing);
      if (model == nullptr) {
         return nullptr;
      }
      uno_set_variables_lower_bounds(model, variables_lower_bounds);
      uno_set_variables_upper_bounds(model, variables_upper_bounds);
      return model;
   });
}

void* uno_create_unconstrained_model(const char* problem_type, uno_int number_variables, uno_int base_indexing) UNO_NOEXCEPT {
   return guard(static_cast<void*>(nullptr), [&]() -> void* {
      if (number_variables <= 0) {
         WARNING << "Please specify a positive number of variables."  << std::endl;
         return nullptr;
      }
      if (problem_type == nullptr || (strcmp(problem_type, UNO_PROBLEM_LINEAR) != 0 &&
            strcmp(problem_type, UNO_PROBLEM_QUADRATIC) != 0 && strcmp(problem_type, UNO_PROBLEM_NONLINEAR) != 0)) {
         WARNING << "Please specify a valid problem type."  << std::endl;
         return nullptr;
      }
      if (base_indexing != UNO_ZERO_BASED_INDEXING && base_indexing != UNO_ONE_BASED_INDEXING) {
         WARNING << "Please specify a valid base indexing." << std::endl;
         return nullptr;
      }
      // owned by a unique_ptr until the model is complete
      auto user_model = std::make_unique<CUserModel>(problem_type, number_variables, base_indexing);
      // set (-inf, +inf) bounds
      const size_t unsigned_number_variables = static_cast<size_t>(number_variables);
      user_model->variables_lower_bounds.resize(unsigned_number_variables, -Inf);
      user_model->variables_upper_bounds.resize(unsigned_number_variables, Inf);
      // create the initial primal point
      user_model->initial_primal_iterate.resize(unsigned_number_variables);
      for (size_t variable_index: Range(unsigned_number_variables)) {
         user_model->initial_primal_iterate[variable_index] = 0.;
      }
      return user_model.release();
   });
}

bool uno_set_model_name(void* model, const char* name) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      user_model->name = (name != nullptr) ? name : "";
      return true;
   });
}

bool uno_set_variables_lower_bounds(void* model, const double* variables_lower_bounds) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // copy the bounds internally
      if (variables_lower_bounds != nullptr) {
         for (size_t variable_index: Range(static_cast<size_t>(user_model->number_variables))) {
            user_model->variables_lower_bounds[variable_index] = variables_lower_bounds[variable_index];
         }
      }
      return true;
   });
}

bool uno_set_variables_upper_bounds(void* model, const double* variables_upper_bounds) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // copy the bounds internally
      if (variables_upper_bounds != nullptr) {
         for (size_t variable_index: Range(static_cast<size_t>(user_model->number_variables))) {
            user_model->variables_upper_bounds[variable_index] = variables_upper_bounds[variable_index];
         }
      }

      return true;
   });
}

bool uno_set_variable_lower_bound(void* model, uno_int variable_index, double lower_bound) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // shift the index with the base indexing
      variable_index -= user_model->base_indexing;
      if (variable_index < 0 || variable_index >= user_model->number_variables) {
         WARNING << "Please specify a valid index."  << std::endl;
         return false;
      }
      user_model->variables_lower_bounds[static_cast<size_t>(variable_index)] = lower_bound;
      return true;
   });
}

bool uno_set_variable_upper_bound(void* model, uno_int variable_index, double upper_bound) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // shift the index with the base indexing
      variable_index -= user_model->base_indexing;
      if (variable_index < 0 || variable_index >= user_model->number_variables) {
         WARNING << "Please specify a valid index."  << std::endl;
         return false;
      }
      user_model->variables_upper_bounds[static_cast<size_t>(variable_index)] = upper_bound;
      return true;
   });
}

bool uno_set_objective(void* model, uno_int optimization_sense, uno_objective_callback objective_function,
      uno_objective_gradient_callback objective_gradient) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
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
   });
}

bool uno_set_constraints(void* model, uno_int number_constraints, uno_constraints_callback constraint_functions,
      const double* constraints_lower_bounds, const double* constraints_upper_bounds, uno_int number_jacobian_nonzeros,
      const uno_int* jacobian_row_indices, const uno_int* jacobian_column_indices, uno_constraints_jacobian_callback jacobian) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (number_constraints <= 0) {
         WARNING << "Please specify a positive number of constraints."  << std::endl;
         return false;
      }
      if (number_jacobian_nonzeros < 0) {
         WARNING << "Please specify a nonnegative number of Jacobian nonzeros."  << std::endl;
         return false;
      }
      if (0 < number_jacobian_nonzeros && (jacobian_row_indices == nullptr || jacobian_column_indices == nullptr)) {
         WARNING << "Please specify valid Jacobian row and column indices."  << std::endl;
         return false;
      }
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // build the new data first: if an allocation fails, the model is left untouched
      const size_t unsigned_number_constraints = static_cast<size_t>(number_constraints);
      const size_t unsigned_number_jacobian_nonzeros = static_cast<size_t>(number_jacobian_nonzeros);
      // copy the bounds internally
      std::vector<double> lower_bounds = user_model->constraints_lower_bounds;
      std::vector<double> upper_bounds = user_model->constraints_upper_bounds;
      lower_bounds.resize(unsigned_number_constraints);
      upper_bounds.resize(unsigned_number_constraints);
      for (size_t constraint_index: Range(unsigned_number_constraints)) {
         if (constraints_lower_bounds != nullptr) {
            lower_bounds[constraint_index] = constraints_lower_bounds[constraint_index];
         }
         if (constraints_upper_bounds != nullptr) {
            upper_bounds[constraint_index] = constraints_upper_bounds[constraint_index];
         }
      }
      // copy the Jacobian sparsity (in C indexing) to allow the calling code to dispose of its vectors
      std::vector<uno_int> row_indices(unsigned_number_jacobian_nonzeros);
      std::vector<uno_int> column_indices(unsigned_number_jacobian_nonzeros);
      for (size_t index: Range(unsigned_number_jacobian_nonzeros)) {
         row_indices[index] = jacobian_row_indices[index] - user_model->base_indexing;
         column_indices[index] = jacobian_column_indices[index] - user_model->base_indexing;
      }
      // the initial dual point
      std::vector<double> initial_dual_iterate(unsigned_number_constraints, 0.);

      // commit the new data (cannot throw). The counts are set last
      user_model->constraint_functions = constraint_functions;
      user_model->constraints_lower_bounds = std::move(lower_bounds);
      user_model->constraints_upper_bounds = std::move(upper_bounds);
      user_model->jacobian_row_indices = std::move(row_indices);
      user_model->jacobian_column_indices = std::move(column_indices);
      user_model->jacobian = jacobian;
      user_model->initial_dual_iterate = std::move(initial_dual_iterate);
      user_model->number_jacobian_nonzeros = number_jacobian_nonzeros;
      user_model->number_constraints = number_constraints;
      return true;
   });
}

bool uno_set_constraints_lower_bounds(void* model, const double* constraints_lower_bounds) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // copy the bounds internally
      if (constraints_lower_bounds != nullptr) {
         for (size_t constraint_index: Range(static_cast<size_t>(user_model->number_constraints))) {
            user_model->constraints_lower_bounds[constraint_index] = constraints_lower_bounds[constraint_index];
         }
      }
      return true;
   });
}

bool uno_set_constraints_upper_bounds(void* model, const double* constraints_upper_bounds) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // copy the bounds internally
      if (constraints_upper_bounds != nullptr) {
         for (size_t constraint_index: Range(static_cast<size_t>(user_model->number_constraints))) {
            user_model->constraints_upper_bounds[constraint_index] = constraints_upper_bounds[constraint_index];
         }
      }
      return true;
   });
}

bool uno_set_constraint_lower_bound(void* model, uno_int constraint_index, double lower_bound) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // shift the index with the base indexing
      constraint_index -= user_model->base_indexing;
      if (constraint_index < 0 || constraint_index >= user_model->number_constraints) {
         WARNING << "Please specify a valid index."  << std::endl;
         return false;
      }
      user_model->constraints_lower_bounds[static_cast<size_t>(constraint_index)] = lower_bound;
      return true;
   });
}

bool uno_set_constraint_upper_bound(void* model, uno_int constraint_index, double upper_bound) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // shift the index with the base indexing
      constraint_index -= user_model->base_indexing;
      if (constraint_index < 0 || constraint_index >= user_model->number_constraints) {
         WARNING << "Please specify a valid index."  << std::endl;
         return false;
      }
      user_model->constraints_upper_bounds[static_cast<size_t>(constraint_index)] = upper_bound;
      return true;
   });
}

bool uno_set_jacobian_operator(void* model, uno_constraints_jacobian_operator_callback jacobian_operator) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      user_model->jacobian_operator = jacobian_operator;
      return true;
   });
}

bool uno_set_jacobian_transposed_operator(void* model, uno_constraints_jacobian_transposed_operator_callback jacobian_transposed_operator) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      user_model->jacobian_transposed_operator = jacobian_transposed_operator;
      return true;
   });
}

bool uno_set_lagrangian_hessian(void* model, uno_int number_hessian_nonzeros, char hessian_triangular_part,
      const uno_int* hessian_row_indices, const uno_int* hessian_column_indices, uno_lagrangian_hessian_callback lagrangian_hessian) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (number_hessian_nonzeros < 0) {
         WARNING << "Please specify a nonnegative number of Lagrangian Hessian nonzeros."  << std::endl;
         return false;
      }
      if (0 < number_hessian_nonzeros && (hessian_row_indices == nullptr || hessian_column_indices == nullptr)) {
         WARNING << "Please specify valid Lagrangian Hessian row and column indices."  << std::endl;
         return false;
      }
      if (hessian_triangular_part != UNO_LOWER_TRIANGLE && hessian_triangular_part != UNO_UPPER_TRIANGLE) {
         WARNING << "Please specify a correct Hessian triangle in {'" << UNO_LOWER_TRIANGLE << "', '" <<
            UNO_UPPER_TRIANGLE << "'}."  << std::endl;
         return false;
      }
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // if the Hessian is empty, the problem is an LP
      if (number_hessian_nonzeros == 0) {
         user_model->problem_type = ProblemType::LINEAR;
      }
      else {
         // copy the Hessian sparsity to allow the calling code to dispose of its vectors. They are built first: if an
         // allocation fails, the model is left untouched
         // from now on, we only maintain the lower triangle
         const size_t unsigned_number_hessian_nonzeros = static_cast<size_t>(number_hessian_nonzeros);
         std::vector<uno_int> row_indices(unsigned_number_hessian_nonzeros);
         std::vector<uno_int> column_indices(unsigned_number_hessian_nonzeros);
         const bool lower_triangle = (hessian_triangular_part == UNO_LOWER_TRIANGLE);
         for (size_t index: Range(unsigned_number_hessian_nonzeros)) {
            row_indices[index] = lower_triangle ? hessian_row_indices[index] : hessian_column_indices[index];
            column_indices[index] = lower_triangle ? hessian_column_indices[index] : hessian_row_indices[index];
         }
         // commit the new data (cannot throw)
         user_model->hessian_row_indices = std::move(row_indices);
         user_model->hessian_column_indices = std::move(column_indices);
         user_model->hessian_triangular_part = UNO_LOWER_TRIANGLE;
         user_model->lagrangian_hessian = lagrangian_hessian;
      }
      // the count is set last
      user_model->number_hessian_nonzeros = number_hessian_nonzeros;
      return true;
   });
}

bool uno_set_lagrangian_hessian_operator(void* model, uno_lagrangian_hessian_operator_callback lagrangian_hessian_operator) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      user_model->lagrangian_hessian_operator = lagrangian_hessian_operator;
      return true;
   });
}

bool uno_set_lagrangian_sign_convention(void* model, uno_int lagrangian_sign_convention) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
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
      user_model->lagrangian_sign_convention = lagrangian_sign_convention;
      return true;
   });
}

bool uno_set_user_data(void* model, void* user_data) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      user_model->user_data = user_data;
      return true;
   });
}

bool uno_set_initial_primal_iterate_component(void* model, uno_int index, double initial_primal_component) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // shift the index with the base indexing
      index -= user_model->base_indexing;
      if (0 <= index && index < user_model->number_variables) {
         user_model->initial_primal_iterate[static_cast<size_t>(index)] = initial_primal_component;
         return true;
      }
      return false;
   });
}

bool uno_set_initial_dual_iterate_component(void* model, uno_int index, double initial_dual_component) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (model == nullptr) {
         WARNING << "Please specify a valid model."  << std::endl;
         return false;
      }
      CUserModel* user_model = static_cast<CUserModel*>(model);
      // shift the index with the base indexing
      index -= user_model->base_indexing;
      if (0 <= index && index < user_model->number_constraints) {
         user_model->initial_dual_iterate[static_cast<size_t>(index)] = initial_dual_component;
         return true;
      }
      return false;
   });
}

bool uno_set_initial_primal_iterate(void* model, const double* initial_primal_iterate) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
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
   });
}

bool uno_set_initial_dual_iterate(void* model, const double* initial_dual_iterate) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
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
   });
}

void* uno_create_solver() UNO_NOEXCEPT {
   return guard(static_cast<void*>(nullptr), [&]() -> void* {
      auto solver = std::make_unique<Solver>();
      DefaultOptions::load(solver->default_options);
      return solver.release();
   });
}

namespace {
   // runs a function that deals with options and converts its exceptions into a status code: no exception crosses the
   // C interface
   template <typename Function>
   uno_int guard_options(Function&& function) noexcept {
      try {
         function();
         return UNO_OPTIONS_SUCCESS;
      }
      catch (const OptionNotFound& exception) {
         report_error(exception.what());
         return UNO_OPTIONS_NOT_FOUND;
      }
      catch (const OptionWrongType& exception) {
         report_error(exception.what());
         return UNO_OPTIONS_WRONG_TYPE;
      }
      catch (const InvalidOptionValue& exception) {
         report_error(exception.what());
         return UNO_OPTIONS_INVALID_VALUE;
      }
      catch (const InvalidOptionFile& exception) {
         report_error(exception.what());
         return UNO_OPTIONS_INVALID_FILE;
      }
      catch (const std::exception& exception) {
         report_error(exception.what());
         return UNO_OPTIONS_ERROR;
      }
      catch (...) {
         report_error("An unknown error occurred");
         return UNO_OPTIONS_ERROR;
      }
   }

   // the options set by the preset (none for "auto", whose options depend on the model)
   Options compute_preset_options(const Options& user_options, const Options& default_options) {
      Options preset_options{};
      const std::string preset = user_options.get_string_optional("preset").value_or(default_options.get_string("preset"));
      if (preset != "auto") {
         Presets::set(preset_options, preset);
      }
      return preset_options;
   }

   // modifies a copy of the user options, then commits it: if any step fails, the options of the solver are unchanged
   template <typename Modifier>
   uno_int modify_user_options(void* solver, Modifier&& modify) noexcept {
      if (solver == nullptr) {
         report_error("Please specify a valid solver.");
         return UNO_OPTIONS_INVALID_ARGUMENT;
      }
      Solver* uno_solver = static_cast<Solver*>(solver);
      return guard_options([&]() {
         Options user_options = uno_solver->user_options;
         modify(user_options);
         Options preset_options = compute_preset_options(user_options, uno_solver->default_options);
         // commit
         uno_solver->user_options = std::move(user_options);
         uno_solver->preset_options = std::move(preset_options);
         uno_solver->options_consistent_with_solve = false;
         uno_solver->solve_options = Options{};
      });
   }

   // sets a particular option
   template <typename Setter>
   uno_int set_user_option(void* solver, const char* option_name, Setter&& setter) noexcept {
      if (option_name == nullptr) {
         report_error("Please specify a valid option name.");
         return UNO_OPTIONS_INVALID_ARGUMENT;
      }
      return modify_user_options(solver, [&](Options& options) {
         setter(options, std::string{option_name});
      });
   }
} // namespace

uno_int uno_set_solver_integer_option(void* solver, const char* option_name, uno_int option_value) UNO_NOEXCEPT {
   return set_user_option(solver, option_name, [&](Options& options, const std::string& name) {
      options.set_integer(name, option_value);
   });
}

uno_int uno_set_solver_double_option(void* solver, const char* option_name, double option_value) UNO_NOEXCEPT {
   return set_user_option(solver, option_name, [&](Options& options, const std::string& name) {
      options.set_double(name, option_value);
   });
}

uno_int uno_set_solver_bool_option(void* solver, const char* option_name, bool option_value) UNO_NOEXCEPT {
   return set_user_option(solver, option_name, [&](Options& options, const std::string& name) {
      options.set_bool(name, option_value);
   });
}

uno_int uno_set_solver_string_option(void* solver, const char* option_name, const char* option_value) UNO_NOEXCEPT {
   if (option_value == nullptr) {
      report_error("Please specify a valid option value.");
      return UNO_OPTIONS_INVALID_ARGUMENT;
   }
   return set_user_option(solver, option_name, [&](Options& options, const std::string& name) {
      options.set_string(name, std::string{option_value});
   });
}

uno_int uno_get_solver_option_type(void* solver, const char* option_name) UNO_NOEXCEPT {
   if (solver == nullptr || option_name == nullptr) {
      return UNO_OPTION_TYPE_NOT_FOUND;
   }
   try {
      const auto type = Options::option_types.find(option_name);
      if (type == Options::option_types.end()) {
         return UNO_OPTION_TYPE_NOT_FOUND;
      }
      return static_cast<uno_int>(type->second);
   }
   catch (...) {
      return UNO_OPTION_TYPE_NOT_FOUND;
   }
}

uno_options_iterator uno_options_begin_iterator() UNO_NOEXCEPT {
   if (Options::option_types.begin() == Options::option_types.end()) {
      // empty list
      return nullptr;
   }
   try {
      return new OptionIterator(Options::option_types.begin());
   }
   catch (...) {
      return nullptr;
   }
}

uno_options_iterator uno_options_end_iterator() UNO_NOEXCEPT {
   // nullptr is the sentinel for end iterator
   return nullptr;
}

void uno_options_next_iterator(uno_options_iterator* iterator) UNO_NOEXCEPT {
   if (iterator == nullptr) {
      report_error("Please specify a valid options iterator.");
      return;
   }
   if (*iterator == nullptr) {
      // do not increment if end has been reached
      return;
   }
   auto& uno_iterator = *static_cast<OptionIterator*>(*iterator);
   ++uno_iterator;
   if (uno_iterator == Options::option_types.end()) {
      // end reached: destroy
      uno_options_destroy_iterator(iterator);
   }
}

const char* uno_options_iterator_name(uno_options_iterator iterator) UNO_NOEXCEPT {
   if (iterator == nullptr) {
      // nullptr if iterator is null
      return nullptr;
   }
   const auto& uno_iterator = *static_cast<OptionIterator*>(iterator);
   return uno_iterator->first.c_str();
}

uno_int uno_options_iterator_type(uno_options_iterator iterator) UNO_NOEXCEPT {
   if (iterator == nullptr) {
      // option not found if iterator is null
      return UNO_OPTION_TYPE_NOT_FOUND;
   }
   const auto& uno_iterator = *static_cast<OptionIterator*>(iterator);
   return static_cast<uno_int>(uno_iterator->second);
}

void uno_options_destroy_iterator(uno_options_iterator* iterator) UNO_NOEXCEPT {
   if (iterator == nullptr || *iterator == nullptr) {
      return;
   }
   delete static_cast<OptionIterator*>(*iterator);
   *iterator = nullptr;
}

uno_int uno_load_solver_option_file(void* solver, const char* file_name) UNO_NOEXCEPT {
   if (file_name == nullptr) {
      report_error("Please specify a valid option file name.");
      return UNO_OPTIONS_INVALID_ARGUMENT;
   }
   return modify_user_options(solver, [&](Options& options) {
      Options::load_option_file(options, file_name);
   });
}

uno_int uno_set_solver_preset(void* solver, const char* preset_name) UNO_NOEXCEPT {
   if (preset_name == nullptr) {
      report_error("Please specify a valid preset.");
      return UNO_OPTIONS_INVALID_ARGUMENT;
   }
   return modify_user_options(solver, [&](Options& options) {
      options.set_string("preset", std::string{preset_name});
   });
}

bool uno_set_solver_callbacks(void* solver, uno_notify_acceptable_iterate_callback notify_acceptable_iterate_callback,
      uno_termination_callback termination_callback, void* user_data) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      if (solver == nullptr) {
         WARNING << "Please specify a valid solver."  << std::endl;
         return false;
      }
      Solver* uno_solver = static_cast<Solver*>(solver);
      // replaces (and deletes) the previous callbacks
      uno_solver->user_callbacks = std::make_unique<CUserCallbacks>(notify_acceptable_iterate_callback,
         termination_callback, user_data);
      return true;
   });
}

bool uno_set_logger_stream_callback(uno_logger_stream_callback logger_stream_callback, void* user_data) UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      // the logger never refers to a deleted stream
      COStream* stream = new COStream(logger_stream_callback, user_data);
      Logger::set_stream(*stream);
      delete c_ostream;
      c_ostream = stream;
      return true;
   });
}

bool uno_reset_logger_stream() UNO_NOEXCEPT {
   return guard(false, [&]() -> bool {
      Logger::set_stream(std::cout);
      delete c_ostream;
      c_ostream = nullptr;
      return true;
   });
}

uno_int uno_optimize(void* solver, void* model) UNO_NOEXCEPT {
   if (solver == nullptr) {
      report_error("Please specify a valid solver.");
      return UNO_OPTIMIZE_INVALID_ARGUMENT;
   }
   Solver* uno_solver = static_cast<Solver*>(solver);
   // the result of a previous optimization does not correspond to this optimization
   uno_solver->result.reset();
   uno_solver->options_consistent_with_solve = false;
   uno_solver->solve_options = Options{};

   if (model == nullptr) {
      report_error("Please specify a valid model.");
      return UNO_OPTIMIZE_INVALID_ARGUMENT;
   }
   CUserModel* user_model = static_cast<CUserModel*>(model);
   if (!user_model->objective_function && !user_model->constraint_functions) {
      report_error("Please specify at least an objective or constraints.");
      return UNO_OPTIMIZE_INVALID_MODEL;
   }

   // the status returned if an exception is raised: it depends on the current stage
   uno_int failure_status = UNO_OPTIMIZE_ERROR;
   try {
      // create an instance of UnoModel, a subclass of Model
      const UnoModel uno_model(*user_model);

      // defaults -> preset -> user options
      failure_status = UNO_OPTIMIZE_INVALID_OPTIONS;
      Options full_options = uno_solver->default_options;
      const std::string preset = uno_solver->user_options.get_string_optional("preset").value_or(uno_solver->default_options.get_string("preset"));
      if (preset == "auto") {
         Presets::set(uno_model, full_options, preset);
      }
      else {
         full_options.overwrite(uno_solver->preset_options);
      }
      full_options.overwrite(uno_solver->user_options);

      // resolve logger
      Logger::set_logger(full_options.get_string("logger"));

      uno_solver->user_options.print("User options");

      // solve the model. Uno::solve reports its errors in the result
      failure_status = UNO_OPTIMIZE_ERROR;
      Result result = uno_solver->solver->solve(uno_model, full_options, *uno_solver->user_callbacks);
      // move the new result into uno_solver
      uno_solver->result = std::make_unique<Result>(std::move(result));
      // keep the options used by the solve, so that the getters can return them (including the values picked by the
      // "auto" preset)
      uno_solver->solve_options = std::move(full_options);
      uno_solver->options_consistent_with_solve = true;
      // flush the logger
      Logger::flush();
      return UNO_OPTIMIZE_SUCCESS;
   }
   catch (const std::bad_alloc& exception) {
      report_error(exception.what());
      failure_status = UNO_OPTIMIZE_ERROR;
   }
   catch (const std::exception& exception) {
      report_error(exception.what());
   }
   catch (...) {
      report_error("An unknown error occurred");
      failure_status = UNO_OPTIMIZE_ERROR;
   }
   try {
      Logger::flush();
   }
   catch (...) {
   }
   return failure_status;
}

namespace {
   // gets the value of an option of type "type": the value set by the user, otherwise by the preset, otherwise the default
   // value. "read" reads the value from a set of options that contains it
   template <typename T, typename Reader>
   uno_int get_solver_option(void* solver, const char* option_name, T* option_value, OptionType type, Reader&& read) noexcept {
      if (solver == nullptr) {
         report_error("Please specify a valid solver.");
         return UNO_OPTIONS_INVALID_ARGUMENT;
      }
      if (option_name == nullptr) {
         report_error("Please specify a valid option name.");
         return UNO_OPTIONS_INVALID_ARGUMENT;
      }
      if (option_value == nullptr) {
         report_error("Please specify a valid pointer to the option value.");
         return UNO_OPTIONS_INVALID_ARGUMENT;
      }
      const Solver* uno_solver = static_cast<const Solver*>(solver);
      uno_int status = UNO_OPTIONS_NOT_SET;
      const uno_int guard_status = guard_options([&]() {
         const std::string name{option_name};
         // throws OptionNotFound if the option does not exist
         if (Options::get_option_type(name) != type) {
            throw OptionWrongType("The option with name " + name + " has a different type");
         }
         // after a solve (and until the options are modified), the options used by the solve. Otherwise, the user
         // options, then the preset options, then the default options
         const std::array<const Options*, 3> sources = uno_solver->options_consistent_with_solve ?
            std::array<const Options*, 3>{&uno_solver->solve_options, nullptr, nullptr} :
            std::array<const Options*, 3>{&uno_solver->user_options, &uno_solver->preset_options, &uno_solver->default_options};
         for (const Options* options: sources) {
            if (options != nullptr && options->contains(name)) {
               *option_value = read(*options, name);
               status = UNO_OPTIONS_SUCCESS;
               return;
            }
         }
      });
      return (guard_status == UNO_OPTIONS_SUCCESS) ? status : guard_status;
   }
} // namespace

uno_int uno_get_solver_integer_option(void* solver, const char* option_name, uno_int* option_value) UNO_NOEXCEPT {
   return get_solver_option(solver, option_name, option_value, OptionType::INTEGER,
      [](const Options& options, const std::string& name) { return options.get_int(name); });
}

uno_int uno_get_solver_double_option(void* solver, const char* option_name, double* option_value) UNO_NOEXCEPT {
   return get_solver_option(solver, option_name, option_value, OptionType::DOUBLE,
      [](const Options& options, const std::string& name) { return options.get_double(name); });
}

uno_int uno_get_solver_bool_option(void* solver, const char* option_name, bool* option_value) UNO_NOEXCEPT {
   return get_solver_option(solver, option_name, option_value, OptionType::BOOL,
      [](const Options& options, const std::string& name) { return options.get_bool(name); });
}

namespace {
   // after a solve, the value of a string option overridden by the solver, otherwise nullptr
   const char* find_option_override(const Solver& solver, const char* option_name) {
      if (solver.result != nullptr && solver.options_consistent_with_solve) {
         for (const auto& option_override: solver.result->option_overrides) {
            if (option_override.option_name == option_name) {
               return option_override.new_value.c_str();
            }
         }
      }
      return nullptr;
   }
} // namespace

uno_int uno_get_solver_string_option(void* solver, const char* option_name, const char** option_value) UNO_NOEXCEPT {
   // the pointer to the stored string (not a copy)
   const uno_int status = get_solver_option(solver, option_name, option_value, OptionType::STRING,
      [](const Options& options, const std::string& name) { return options.get_string(name).c_str(); });
   // after a solve, an option overridden by the solver has the value used by the solver. The arguments are valid if
   // the status is UNO_OPTIONS_SUCCESS or UNO_OPTIONS_NOT_SET
   if (status == UNO_OPTIONS_SUCCESS || status == UNO_OPTIONS_NOT_SET) {
      if (const char* overridden_value = find_option_override(*static_cast<const Solver*>(solver), option_name)) {
         *option_value = overridden_value;
         return UNO_OPTIONS_SUCCESS;
      }
   }
   return status;
}

namespace {
   // checks the arguments and gets the result of the last optimization. If there is no result (e.g. before the first
   // optimization), returns UNO_RESULT_NOT_AVAILABLE without a message: this is a state, not an error
   uno_int get_result(void* solver, bool valid_output, const Result*& result) noexcept {
      if (solver == nullptr) {
         report_error("Please specify a valid solver.");
         return UNO_RESULT_INVALID_ARGUMENT;
      }
      if (!valid_output) {
         report_error("Please specify a valid pointer to the value.");
         return UNO_RESULT_INVALID_ARGUMENT;
      }
      result = static_cast<const Solver*>(solver)->result.get();
      return (result != nullptr) ? UNO_RESULT_SUCCESS : UNO_RESULT_NOT_AVAILABLE;
   }

   // gets a member of the result of the last optimization
   template <typename T, typename Member>
   uno_int get_result_value(void* solver, T* value, const Member Result::* member) noexcept {
      const Result* result{nullptr};
      const uno_int status = get_result(solver, value != nullptr, result);
      if (status == UNO_RESULT_SUCCESS) {
         *value = static_cast<T>(result->*member);
      }
      return status;
   }

   // gets one component of a vector of the result, whose size is "size"
   uno_int get_result_component(void* solver, uno_int index, double* value, Vector<double> Result::* vector,
         const size_t Result::* size) noexcept {
      const Result* result{nullptr};
      const uno_int status = get_result(solver, value != nullptr, result);
      if (status != UNO_RESULT_SUCCESS) {
         return status;
      }
      // shift the index with the base indexing
      index -= result->base_indexing;
      if (index < 0 || result->*size <= static_cast<size_t>(index)) {
         report_error("The index is not valid.");
         return UNO_RESULT_INVALID_INDEX;
      }
      *value = (result->*vector)[static_cast<size_t>(index)];
      return UNO_RESULT_SUCCESS;
   }

   // copies a vector of the result, whose size is "size". A null pointer is allowed for an empty vector
   uno_int get_result_vector(void* solver, double* values, Vector<double> Result::* vector,
         const size_t Result::* size) noexcept {
      const Result* result{nullptr};
      const uno_int status = get_result(solver, true, result);
      if (status != UNO_RESULT_SUCCESS) {
         return status;
      }
      const size_t number_values = result->*size;
      if (values == nullptr && 0 < number_values) {
         report_error("Please specify a valid pointer to the values.");
         return UNO_RESULT_INVALID_ARGUMENT;
      }
      const Vector<double>& result_vector = result->*vector;
      for (size_t index: Range(number_values)) {
         values[index] = result_vector[index];
      }
      return UNO_RESULT_SUCCESS;
   }
} // namespace

uno_int uno_get_method_description(void* solver, const char** method_description) UNO_NOEXCEPT {
   const Result* result{nullptr};
   const uno_int status = get_result(solver, method_description != nullptr, result);
   if (status == UNO_RESULT_SUCCESS) {
      *method_description = static_cast<const Solver*>(solver)->solver->get_method_description().c_str();
   }
   return status;
}

uno_int uno_get_optimization_status(void* solver, uno_int* optimization_status) UNO_NOEXCEPT {
   return get_result_value(solver, optimization_status, &Result::optimization_status);
}

uno_int uno_get_solution_status(void* solver, uno_int* solution_status) UNO_NOEXCEPT {
   return get_result_value(solver, solution_status, &Result::solution_status);
}

uno_int uno_get_solution_objective(void* solver, double* solution_objective) UNO_NOEXCEPT {
   return get_result_value(solver, solution_objective, &Result::solution_objective);
}

uno_int uno_get_primal_solution_component(void* solver, uno_int index, double* primal_component) UNO_NOEXCEPT {
   return get_result_component(solver, index, primal_component, &Result::primal_solution, &Result::number_variables);
}

uno_int uno_get_constraint_dual_solution_component(void* solver, uno_int index, double* constraint_dual_component) UNO_NOEXCEPT {
   return get_result_component(solver, index, constraint_dual_component, &Result::constraint_dual_solution,
      &Result::number_constraints);
}

uno_int uno_get_lower_bound_dual_solution_component(void* solver, uno_int index, double* lower_bound_dual_component) UNO_NOEXCEPT {
   return get_result_component(solver, index, lower_bound_dual_component, &Result::lower_bound_dual_solution,
      &Result::number_variables);
}

uno_int uno_get_upper_bound_dual_solution_component(void* solver, uno_int index, double* upper_bound_dual_component) UNO_NOEXCEPT {
   return get_result_component(solver, index, upper_bound_dual_component, &Result::upper_bound_dual_solution,
      &Result::number_variables);
}

uno_int uno_get_primal_solution(void* solver, double* primal_solution) UNO_NOEXCEPT {
   return get_result_vector(solver, primal_solution, &Result::primal_solution, &Result::number_variables);
}

uno_int uno_get_constraint_dual_solution(void* solver, double* constraint_dual_solution) UNO_NOEXCEPT {
   return get_result_vector(solver, constraint_dual_solution, &Result::constraint_dual_solution,
      &Result::number_constraints);
}

uno_int uno_get_lower_bound_dual_solution(void* solver, double* lower_bound_dual_solution) UNO_NOEXCEPT {
   return get_result_vector(solver, lower_bound_dual_solution, &Result::lower_bound_dual_solution,
      &Result::number_variables);
}

uno_int uno_get_upper_bound_dual_solution(void* solver, double* upper_bound_dual_solution) UNO_NOEXCEPT {
   return get_result_vector(solver, upper_bound_dual_solution, &Result::upper_bound_dual_solution,
      &Result::number_variables);
}

uno_int uno_get_solution_constraints(void* solver, double* constraint_values) UNO_NOEXCEPT {
   return get_result_vector(solver, constraint_values, &Result::constraint_values, &Result::number_constraints);
}

uno_int uno_get_solution_primal_feasibility(void* solver, double* primal_feasibility) UNO_NOEXCEPT {
   return get_result_value(solver, primal_feasibility, &Result::solution_primal_feasibility);
}

uno_int uno_get_solution_stationarity(void* solver, double* stationarity) UNO_NOEXCEPT {
   return get_result_value(solver, stationarity, &Result::solution_stationarity);
}

uno_int uno_get_solution_complementarity(void* solver, double* complementarity) UNO_NOEXCEPT {
   return get_result_value(solver, complementarity, &Result::solution_complementarity);
}

uno_int uno_get_number_iterations(void* solver, uno_int* number_iterations) UNO_NOEXCEPT {
   return get_result_value(solver, number_iterations, &Result::number_iterations);
}

uno_int uno_get_cpu_time(void* solver, double* cpu_time) UNO_NOEXCEPT {
   const Result* result{nullptr};
   const uno_int status = get_result(solver, cpu_time != nullptr, result);
   if (status == UNO_RESULT_SUCCESS) {
      *cpu_time = result->timers.solve.get_duration();
   }
   return status;
}

uno_int uno_get_number_objective_evaluations(void* solver, uno_int* number_evaluations) UNO_NOEXCEPT {
   return get_result_value(solver, number_evaluations, &Result::number_objective_evaluations);
}

uno_int uno_get_number_constraint_evaluations(void* solver, uno_int* number_evaluations) UNO_NOEXCEPT {
   return get_result_value(solver, number_evaluations, &Result::number_constraint_evaluations);
}

uno_int uno_get_number_objective_gradient_evaluations(void* solver, uno_int* number_evaluations) UNO_NOEXCEPT {
   return get_result_value(solver, number_evaluations, &Result::number_objective_gradient_evaluations);
}

uno_int uno_get_number_jacobian_evaluations(void* solver, uno_int* number_evaluations) UNO_NOEXCEPT {
   return get_result_value(solver, number_evaluations, &Result::number_jacobian_evaluations);
}

uno_int uno_get_number_hessian_evaluations(void* solver, uno_int* number_evaluations) UNO_NOEXCEPT {
   return get_result_value(solver, number_evaluations, &Result::number_hessian_evaluations);
}

uno_int uno_get_number_subproblem_solved_evaluations(void* solver, uno_int* number_subproblems) UNO_NOEXCEPT {
   return get_result_value(solver, number_subproblems, &Result::number_subproblems_solved);
}

void uno_destroy_model(void* model) UNO_NOEXCEPT {
   if (model != nullptr) {
      delete static_cast<CUserModel*>(model);
   }
}

void uno_destroy_solver(void* solver) UNO_NOEXCEPT {
   delete static_cast<Solver*>(solver);
}
