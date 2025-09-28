// Copyright (c) 2025 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <iostream>
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

using namespace uno;

using CUserModel = UserModel<Objective, ObjectiveGradient, Constraints, Jacobian, JacobianOperator, JacobianTransposedOperator,
   Hessian, HessianOperator, const double*, void*>;

// UnoModel contains an instance of UserModel and complies with the Model interface
class UnoModel: public Model {
public:
   explicit UnoModel(const CUserModel& user_model):
         Model("C model", static_cast<size_t>(user_model.number_variables), static_cast<size_t>(user_model.number_constraints),
            static_cast<double>(user_model.optimization_sense)),
         user_model(user_model),
         equality_constraints_collection(this->equality_constraints),
         inequality_constraints_collection(this->inequality_constraints) {
      this->find_fixed_variables(this->fixed_variables);
      this->partition_constraints(this->equality_constraints, this->inequality_constraints);
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
         const int32_t return_code = this->user_model.objective_function(this->user_model.number_variables, x.data(),
            &objective_value, this->user_model.user_data);
         if (0 < return_code) {
            throw FunctionEvaluationError();
         }
         objective_value *= this->optimization_sense;
      }
      return objective_value;
   }

   void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override {
      if (this->user_model.constraint_functions != nullptr) {
         const int32_t return_code = this->user_model.constraint_functions(this->user_model.number_variables,
            this->user_model.number_constraints, x.data(), constraints.data(), this->user_model.user_data);
         if (0 < return_code) {
            throw FunctionEvaluationError();
         }
      }
   }

   // dense objective gradient
   void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
      if (this->user_model.objective_gradient != nullptr) {
         const int32_t return_code = this->user_model.objective_gradient(this->user_model.number_variables, x.data(),
            gradient.data(), this->user_model.user_data);
         if (0 < return_code) {
            throw GradientEvaluationError();
         }
         for (size_t variable_index: Range(this->number_variables)) {
            gradient[variable_index] *= this->optimization_sense;
         }
      }
   }

   // sparsity patterns of Jacobian and Hessian
   void compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
         MatrixOrder /*matrix_order*/) const override {
      // copy the indices of the user sparsity patterns to the Uno vectors
      std::copy_n(this->user_model.jacobian_row_indices.data(), static_cast<size_t>(this->user_model.number_jacobian_nonzeros), row_indices);
      std::copy_n(this->user_model.jacobian_column_indices.data(), static_cast<size_t>(this->user_model.number_jacobian_nonzeros), column_indices);
      // TODO matrix_order

      // handle the solver indexing
      if (this->user_model.base_indexing != solver_indexing) {
         const int indexing_difference = solver_indexing - this->user_model.base_indexing;
         for (size_t index: Range(static_cast<size_t>(this->user_model.number_jacobian_nonzeros))) {
            row_indices[index] += indexing_difference;
            column_indices[index] += indexing_difference;
         }
      }
   }

   void compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const override {
      // copy the indices of the user sparsity patterns to the Uno vectors
      std::copy_n(this->user_model.hessian_row_indices.data(), static_cast<size_t>(this->user_model.number_hessian_nonzeros), row_indices);
      std::copy_n(this->user_model.hessian_column_indices.data(), static_cast<size_t>(this->user_model.number_hessian_nonzeros), column_indices);

      // handle the solver indexing
      if (this->user_model.base_indexing != solver_indexing) {
         const int indexing_difference = solver_indexing - this->user_model.base_indexing;
         for (size_t index: Range(static_cast<size_t>(this->user_model.number_hessian_nonzeros))) {
            row_indices[index] += indexing_difference;
            column_indices[index] += indexing_difference;
         }
      }
   }

   // numerical evaluations of Jacobian and Hessian
   void evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const override {
      if (this->user_model.constraint_jacobian != nullptr) {
         const int32_t return_code = this->user_model.constraint_jacobian(this->user_model.number_variables,
            this->user_model.number_jacobian_nonzeros, x.data(), jacobian_values, this->user_model.user_data);
         if (0 < return_code) {
            throw GradientEvaluationError();
         }
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
         const int32_t return_code = this->user_model.lagrangian_hessian(this->user_model.number_variables,
            this->user_model.number_constraints, this->user_model.number_hessian_nonzeros, x.data(), objective_multiplier,
            multipliers.data(), hessian_values, this->user_model.user_data);
         // flip the signs of the multipliers back
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            const_cast<Vector<double>&>(multipliers).scale(-1.);
         }
         if (0 < return_code) {
            throw HessianEvaluationError();
         }
      }
      else {
         throw std::runtime_error("evaluate_lagrangian_hessian not implemented");
      }
   }

   void compute_jacobian_vector_product(const double* x, const double* vector, double* result) const override {
      if (this->user_model.jacobian_operator != nullptr) {
         const int32_t return_code = this->user_model.jacobian_operator(this->user_model.number_variables,
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
         const int32_t return_code = this->user_model.jacobian_transposed_operator(this->user_model.number_variables,
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
         const int32_t return_code = this->user_model.lagrangian_hessian_operator(this->user_model.number_variables,
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
      if (this->user_model.variables_lower_bounds != nullptr) {
         return this->user_model.variables_lower_bounds[variable_index];
      }
      return -INF<double>;
   }

   [[nodiscard]] double variable_upper_bound(size_t variable_index) const override {
      if (this->user_model.variables_upper_bounds != nullptr) {
         return this->user_model.variables_upper_bounds[variable_index];
      }
      return INF<double>;
   }

   [[nodiscard]] const SparseVector<size_t>& get_slacks() const override {
      return this->slacks;
   }

   [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override {
      return this->fixed_variables;
   }

   [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override {
      if (this->user_model.constraints_lower_bounds != nullptr) {
         return this->user_model.constraints_lower_bounds[constraint_index];
      }
      return -INF<double>;
   }

   [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override {
      if (this->user_model.constraints_upper_bounds != nullptr) {
         return this->user_model.constraints_upper_bounds[constraint_index];
      }
      return INF<double>;
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

   void initial_primal_point(Vector<double>& x) const override {
      if (this->user_model.initial_primal_iterate != nullptr) {
         std::copy_n(this->user_model.initial_primal_iterate, this->user_model.number_variables, x.data());
      }
      else {
         x.fill(0.);
      }
   }

   void initial_dual_point(Vector<double>& multipliers) const override {
      if (this->user_model.initial_dual_iterate != nullptr) {
         std::copy_n(this->user_model.initial_dual_iterate, this->user_model.number_constraints, multipliers.data());
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            multipliers.scale(-1.);
         }
      }
      else {
         multipliers.fill(0.);
      }
   }

   void postprocess_solution(Iterate& iterate) const override {
      // flip the signs of the multipliers, depending on what the sign convention of the Lagrangian is, and whether
      // we maximize
      iterate.multipliers.constraints *= -this->user_model.lagrangian_sign_convention * this->optimization_sense;
      iterate.multipliers.lower_bounds *= -this->user_model.lagrangian_sign_convention * this->optimization_sense;
      iterate.multipliers.upper_bounds *= -this->user_model.lagrangian_sign_convention * this->optimization_sense;
      iterate.evaluations.objective *= this->optimization_sense;
   }

   [[nodiscard]] size_t number_jacobian_nonzeros() const override {
      return static_cast<size_t>(this->user_model.number_jacobian_nonzeros);
   }

   [[nodiscard]] size_t number_hessian_nonzeros() const override {
      return static_cast<size_t>(this->user_model.number_hessian_nonzeros);
   }

protected:
   const CUserModel& user_model;
   const SparseVector<size_t> slacks{};
   Vector<size_t> fixed_variables{};
   const ForwardRange linear_constraints{0};
   std::vector<size_t> equality_constraints;
   CollectionAdapter<std::vector<size_t>> equality_constraints_collection;
   std::vector<size_t> inequality_constraints;
   CollectionAdapter<std::vector<size_t>> inequality_constraints_collection;
};

struct Solver {
   Uno* solver;
   Options* options;
   Result* result;
};

void uno_get_version(int32_t* major, int32_t* minor, int32_t* patch) {
   *major = UNO_VERSION_MAJOR;
   *minor = UNO_VERSION_MINOR;
   *patch = UNO_VERSION_PATCH;
}

void* uno_create_model(char problem_type, int32_t number_variables, const double* variables_lower_bounds,
      const double* variables_upper_bounds, int32_t base_indexing) {
   if (number_variables <= 0) {
      std::cout << "Please specify a positive number of variables.\n";
      return nullptr;
   }
   if (problem_type != UNO_PROBLEM_LINEAR && problem_type != UNO_PROBLEM_QUADRATIC && problem_type != UNO_PROBLEM_NONLINEAR) {
      std::cout << "Please specify a valid problem type.\n";
      return nullptr;
   }
   if (base_indexing != UNO_ZERO_BASED_INDEXING && base_indexing != UNO_ONE_BASED_INDEXING) {
      std::cout << "Please specify a valid base indexing.\n";
      return nullptr;
   }
   return new CUserModel(problem_type, number_variables, variables_lower_bounds, variables_upper_bounds, base_indexing);
}

bool uno_set_objective(void* model, int32_t optimization_sense, Objective objective_function,
      ObjectiveGradient objective_gradient) {
   if (optimization_sense != UNO_MINIMIZE && optimization_sense != UNO_MAXIMIZE) {
      std::cout << "Please specify a valid objective sense.\n";
      return false;
   }

   assert(model != nullptr);
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->optimization_sense = optimization_sense;
   user_model->objective_function = objective_function;
   user_model->objective_gradient = objective_gradient;
   return true;
}

bool uno_set_constraints(void* model, int32_t number_constraints, Constraints constraint_functions,
      const double* constraints_lower_bounds, const double* constraints_upper_bounds, int32_t number_jacobian_nonzeros,
      const int32_t* jacobian_row_indices, const int32_t* jacobian_column_indices, Jacobian constraint_jacobian) {
   if (number_constraints <= 0) {
      std::cout << "Please specify a positive number of constraints.\n";
      return false;
   }

   assert(model != nullptr);
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->number_constraints = number_constraints;
   user_model->constraint_functions = constraint_functions;
   user_model->constraints_lower_bounds = constraints_lower_bounds;
   user_model->constraints_upper_bounds = constraints_upper_bounds;
   user_model->number_jacobian_nonzeros = number_jacobian_nonzeros;
   // copy the Jacobian sparsity to allow the calling code to dispose of its vectors
   user_model->jacobian_row_indices.resize(static_cast<size_t>(number_jacobian_nonzeros));
   user_model->jacobian_column_indices.resize(static_cast<size_t>(number_jacobian_nonzeros));
   for (size_t index: Range(static_cast<size_t>(number_jacobian_nonzeros))) {
      user_model->jacobian_row_indices[index] = jacobian_row_indices[index];
      user_model->jacobian_column_indices[index] = jacobian_column_indices[index];
   }
   user_model->constraint_jacobian = constraint_jacobian;
   return true;
}

bool uno_set_jacobian_operator(void* model, JacobianOperator jacobian_operator) {
   assert(model != nullptr);
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->jacobian_operator = jacobian_operator;
   return true;
}

bool uno_set_jacobian_transposed_operator(void* model, JacobianTransposedOperator jacobian_transposed_operator) {
   assert(model != nullptr);
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->jacobian_transposed_operator = jacobian_transposed_operator;
   return true;
}

bool uno_set_lagrangian_hessian(void* model, int32_t number_hessian_nonzeros, char hessian_triangular_part,
      const int32_t* hessian_row_indices, const int32_t* hessian_column_indices, Hessian lagrangian_hessian,
      double lagrangian_sign_convention) {
   if (number_hessian_nonzeros <= 0) {
      std::cout << "Please specify a positive number of Lagrangian Hessian nonzeros.\n";
      return false;
   }
   if (hessian_triangular_part != UNO_LOWER_TRIANGLE && hessian_triangular_part != UNO_UPPER_TRIANGLE) {
      std::cout << "Please specify a correct Hessian triangle in {'" << UNO_LOWER_TRIANGLE << "', '" <<
         UNO_UPPER_TRIANGLE << "'}.\n";
      return false;
   }
   if (lagrangian_sign_convention != UNO_MULTIPLIER_NEGATIVE && lagrangian_sign_convention != UNO_MULTIPLIER_POSITIVE) {
      std::cout << "Please specify a correct Lagrangian sign convention in {" << UNO_MULTIPLIER_NEGATIVE << ", " <<
         UNO_MULTIPLIER_POSITIVE << "}.\n";
      return false;
   }

   assert(model != nullptr);
   CUserModel* user_model = static_cast<CUserModel*>(model);
   // make sure that the sign convention is consistent with that of the Hessian operator
   if (user_model->lagrangian_hessian_operator != nullptr && user_model->lagrangian_sign_convention != lagrangian_sign_convention) {
      std::cout << "Please specify a Lagrangian sign convention consistent with that of the Hessian operator.\n";
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

bool uno_set_lagrangian_hessian_operator(void* model, int32_t number_hessian_nonzeros, HessianOperator lagrangian_hessian_operator,
      double lagrangian_sign_convention) {
   if (number_hessian_nonzeros <= 0) {
      std::cout << "Please specify a positive number of Lagrangian Hessian nonzeros.\n";
      return false;
   }
   if (lagrangian_sign_convention != UNO_MULTIPLIER_NEGATIVE && lagrangian_sign_convention != UNO_MULTIPLIER_POSITIVE) {
      std::cout << "Please specify a Lagrangian sign convention in {" << UNO_MULTIPLIER_NEGATIVE << ", " <<
         UNO_MULTIPLIER_POSITIVE << "}.\n";
      return false;
   }

   assert(model != nullptr);
   CUserModel* user_model = static_cast<CUserModel*>(model);
   // make sure that the sign convention is consistent with that of the Hessian function
   if (user_model->lagrangian_hessian != nullptr && user_model->lagrangian_sign_convention != lagrangian_sign_convention) {
      std::cout << "Please specify a Lagrangian sign convention consistent with that of the Hessian function.\n";
      return false;
   }
   user_model->number_hessian_nonzeros = number_hessian_nonzeros;
   user_model->lagrangian_hessian_operator = lagrangian_hessian_operator;
   user_model->lagrangian_sign_convention = lagrangian_sign_convention;
   return true;
}

bool uno_set_user_data(void* model, void* user_data) {
   assert(model != nullptr);
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->user_data = user_data;
   return true;
}

bool uno_set_initial_primal_iterate(void* model, double* initial_primal_iterate) {
   assert(model != nullptr);
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->initial_primal_iterate = initial_primal_iterate;
   return true;
}

bool uno_set_initial_dual_iterate(void* model, double* initial_dual_iterate) {
   assert(model != nullptr);
   CUserModel* user_model = static_cast<CUserModel*>(model);
   user_model->initial_dual_iterate = initial_dual_iterate;
   return true;
}

void uno_set_option(void* options, const char* option_name, const char* option_value) {
   assert(options != nullptr);
   uno::Options* uno_options = static_cast<uno::Options*>(options);
   uno_options->set(option_name, option_value);
}

void* uno_create_solver() {
   // default options
   Options* options = new Options;
   DefaultOptions::load(*options);

   // Uno solver
   Uno* uno_solver = new Uno;
   Solver* solver = new Solver{uno_solver, options, nullptr}; // no result yet
   return solver;
}

void uno_set_solver_option(void* solver, const char* option_name, const char* option_value) {
   Solver* uno_solver = static_cast<Solver*>(solver);
   uno_solver->options->set(option_name, option_value);
}

void uno_set_solver_preset(void* solver, const char* preset_name) {
   Solver* uno_solver = static_cast<Solver*>(solver);
   Presets::set(*uno_solver->options, preset_name);
}

void uno_optimize(void* solver, void* model) {
   // check the model
   assert(model != nullptr);
   CUserModel* user_model = static_cast<CUserModel*>(model);
   if (!user_model->objective_function && !user_model->constraint_functions) {
      std::cout << "Please specify at least an objective or constraints.\n";
      return;
   }

   assert(solver != nullptr);
   Solver* uno_solver = static_cast<Solver*>(solver);

   // create an instance of UnoModel, a subclass of Model, and solve the model using Uno
   const UnoModel uno_model(*user_model);
   Logger::set_logger(uno_solver->options->get_string("logger"));
   Result result = uno_solver->solver->solve(uno_model, *uno_solver->options);
   // clean up the previous result (if any)
   delete uno_solver->result;
   // move the new result into uno_solver
   uno_solver->result = new Result(std::move(result));
}

// auxiliary function
Result* uno_get_result(void* solver) {
   assert(solver != nullptr);
   Solver* uno_solver = static_cast<Solver*>(solver);
   assert(uno_solver->result != nullptr);
   return uno_solver->result;
}

int32_t uno_get_optimization_status(void* solver) {
   Result* result = uno_get_result(solver);
   return static_cast<int32_t>(result->optimization_status);
}

int32_t uno_get_solution_status(void* solver) {
   Result* result = uno_get_result(solver);
   return static_cast<int32_t>(result->solution_status);
}

double uno_get_solution_objective(void* solver) {
   Result* result = uno_get_result(solver);
   return result->solution_objective;
}

void uno_get_primal_solution(void* solver, double* primal_solution) {
   Result* result = uno_get_result(solver);
   std::copy_n(result->primal_solution.data(), result->number_variables, primal_solution);
}

void uno_get_constraint_dual_solution(void* solver, double* constraint_dual_solution) {
   Result* result = uno_get_result(solver);
   std::copy_n(result->constraint_dual_solution.data(), result->number_constraints, constraint_dual_solution);
}

void uno_get_lower_bound_dual_solution(void* solver, double* lower_bound_dual_solution) {
   Result* result = uno_get_result(solver);
   std::copy_n(result->lower_bound_dual_solution.data(), result->number_variables, lower_bound_dual_solution);
}

void uno_get_upper_bound_dual_solution(void* solver, double* upper_bound_dual_solution) {
   Result* result = uno_get_result(solver);
   std::copy_n(result->upper_bound_dual_solution.data(), result->number_variables, upper_bound_dual_solution);
}

double uno_get_solution_primal_feasibility(void* solver) {
   Result* result = uno_get_result(solver);
   return result->solution_primal_feasibility;
}

double uno_get_solution_dual_feasibility(void* solver) {
   Result* result = uno_get_result(solver);
   return result->solution_dual_feasibility;
}

double uno_get_solution_complementarity(void* solver) {
   Result* result = uno_get_result(solver);
   return result->solution_complementarity;
}

void uno_destroy_model(void* model) {
   assert(model != nullptr);
   delete static_cast<CUserModel*>(model);
}

void uno_destroy_solver(void* solver) {
   assert(solver != nullptr);
   Solver* uno_solver = static_cast<Solver*>(solver);
   delete uno_solver->solver;
   delete uno_solver->options;
   delete uno_solver->result;
   delete uno_solver;
}