// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <iostream>
#include "Uno_C_API.h"
#include "Uno.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "model/Model.hpp"
#include "options/DefaultOptions.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

class UserModel {
public:
   UserModel(char problem_type, int32_t number_variables, double* variables_lower_bounds, double* variables_upper_bounds,
      int32_t vector_indexing):
         problem_type(problem_type),
         vector_indexing(vector_indexing),
         number_variables(number_variables),
         variables_lower_bounds(variables_lower_bounds),
         variables_upper_bounds(variables_upper_bounds) {
   }

   ~UserModel() = default;

   const char problem_type; // 'L' for linear, 'Q' for quadratic, 'N' for nonlinear
   const int32_t vector_indexing; // 0 for C-style indexing, 1 for Fortran-style indexing

   // variables
   const int32_t number_variables;
   double* variables_lower_bounds{nullptr};
   double* variables_upper_bounds{nullptr};

   // objective
   double objective_sense{1.};
   Objective objective_function{nullptr};
   ObjectiveGradient objective_gradient{nullptr};

   // constraints
   int32_t number_constraints{0};
   Constraints constraint_functions{nullptr};
   double* constraints_lower_bounds{nullptr};
   double* constraints_upper_bounds{nullptr};
   int32_t number_jacobian_nonzeros{0};
   int32_t* jacobian_row_indices{nullptr};
   int32_t* jacobian_column_indices{nullptr};
   Jacobian constraint_jacobian{nullptr};
   JacobianOperator jacobian_operator{nullptr};
   JacobianTransposedOperator jacobian_transposed_operator{nullptr};

   // Hessian
   int32_t number_hessian_nonzeros{0};
   // lower ('L') or upper ('U')
   char hessian_triangular_part{}; // default is empty
   int32_t* hessian_row_indices{nullptr};
   int32_t* hessian_column_indices{nullptr};
   Hessian lagrangian_hessian{nullptr};
   HessianOperator lagrangian_hessian_operator{nullptr};
   double lagrangian_sign_convention{-1.};

   void* user_data{nullptr};

   // initial iterate
   double* initial_primal_iterate{nullptr};
   double* initial_dual_iterate{nullptr};
};

class UnoModel: public uno::Model {
public:
   explicit UnoModel(const UserModel& user_model):
      Model("C model", static_cast<size_t>(user_model.number_variables), static_cast<size_t>(user_model.number_constraints),
         user_model.objective_sense),
      user_model(user_model) {
   }

   // TODO handle failures

   // function evaluations
   [[nodiscard]] double evaluate_objective(const uno::Vector<double>& x) const override {
      double objective_value{0.};
      if (this->user_model.objective_function != nullptr) {
         this->user_model.objective_function(this->user_model.number_variables, x.data(), &objective_value,
            this->user_model.user_data);
      }
      return objective_value;
   }

   void evaluate_constraints(const uno::Vector<double>& x, std::vector<double>& constraints) const override {
      if (this->user_model.constraint_functions != nullptr) {
         this->user_model.constraint_functions(this->user_model.number_variables, this->user_model.number_constraints,
            x.data(), constraints.data(), this->user_model.user_data);
      }
   }

   // dense objective gradient
   void evaluate_objective_gradient(const uno::Vector<double>& x, uno::Vector<double>& gradient) const override {
      if (this->user_model.objective_gradient != nullptr) {
         this->user_model.objective_gradient(this->user_model.number_variables, x.data(), gradient.data(),
            this->user_model.user_data);
      }
   }

   // sparsity patterns of Jacobian and Hessian
   void compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
      uno::MatrixOrder matrix_format) const override;
   void compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const override;

   // numerical evaluations of Jacobian and Hessian
   void evaluate_constraint_jacobian(const uno::Vector<double>& x, double* jacobian_values) const override;
   void evaluate_lagrangian_hessian(const uno::Vector<double>& x, double objective_multiplier, const uno::Vector<double>& multipliers,
      double* hessian_values) const override;
   void compute_hessian_vector_product(const double* vector, double objective_multiplier, const uno::Vector<double>& multipliers,
      double* result) const override;

   [[nodiscard]] double variable_lower_bound(size_t variable_index) const override {
      if (this->user_model.variables_lower_bounds != nullptr) {
         return this->user_model.variables_lower_bounds[variable_index];
      }
      return -uno::INF<double>;
   }

   [[nodiscard]] double variable_upper_bound(size_t variable_index) const override {
      if (this->user_model.variables_upper_bounds != nullptr) {
         return this->user_model.variables_upper_bounds[variable_index];
      }
      return uno::INF<double>;
   }

   [[nodiscard]] const uno::SparseVector<size_t>& get_slacks() const override {
      return this->slacks;
   }

   [[nodiscard]] const uno::Vector<size_t>& get_fixed_variables() const override;

   [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
   [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;
   [[nodiscard]] const uno::Collection<size_t>& get_equality_constraints() const override;
   [[nodiscard]] const uno::Collection<size_t>& get_inequality_constraints() const override;
   [[nodiscard]] const uno::Collection<size_t>& get_linear_constraints() const override;

   void initial_primal_point(uno::Vector<double>& x) const override;
   void initial_dual_point(uno::Vector<double>& multipliers) const override;
   void postprocess_solution(uno::Iterate& iterate, uno::IterateStatus termination_status) const override;

   [[nodiscard]] size_t number_jacobian_nonzeros() const override;
   [[nodiscard]] size_t number_hessian_nonzeros() const override;

protected:
   const UserModel& user_model;
   const uno::SparseVector<size_t> slacks{0};
};

void uno_get_version(int32_t* major, int32_t* minor, int32_t* patch) {
   *major = uno_version_major;
   *minor = uno_version_minor;
   *patch = uno_version_patch;
}

void* uno_create_model(char problem_type, int32_t number_variables, double* variables_lower_bounds,
      double* variables_upper_bounds, int32_t vector_indexing) {
   if (number_variables <= 0) {
      std::cout << "Please specify a positive number of variables.\n";
      return nullptr;
   }
   if (vector_indexing != 0 && vector_indexing != 1) {
      std::cout << "Please specify a valid vector indexing.\n";
      return nullptr;
   }
   return new UserModel(problem_type, number_variables, variables_lower_bounds, variables_upper_bounds, vector_indexing);
}

void uno_set_objective(void* model, double objective_sense, Objective objective_function, ObjectiveGradient objective_gradient) {
   if (objective_sense != 1. && objective_sense != -1.) {
      std::cout << "Please specify a valid objective sense.\n";
      return;
   }

   assert(model != nullptr);
   UserModel* user_model = static_cast<UserModel*>(model);
   user_model->objective_sense = objective_sense;
   user_model->objective_function = objective_function;
   user_model->objective_gradient = objective_gradient;
}

void uno_set_constraints(void* model, int32_t number_constraints, Constraints constraint_functions,
      double* constraints_lower_bounds, double* constraints_upper_bounds, int32_t number_jacobian_nonzeros,
      int32_t* jacobian_row_indices, int32_t* jacobian_column_indices, Jacobian constraint_jacobian) {
   if (number_constraints <= 0) {
      std::cout << "Please specify a positive number of constraints.\n";
      return;
   }

   assert(model != nullptr);
   UserModel* user_model = static_cast<UserModel*>(model);
   user_model->number_constraints = number_constraints;
   user_model->constraint_functions = constraint_functions;
   user_model->constraints_lower_bounds = constraints_lower_bounds;
   user_model->constraints_upper_bounds = constraints_upper_bounds;
   user_model->number_jacobian_nonzeros = number_jacobian_nonzeros;
   user_model->jacobian_row_indices = jacobian_row_indices;
   user_model->jacobian_column_indices = jacobian_column_indices;
   user_model->constraint_jacobian = constraint_jacobian;
}

void uno_set_jacobian_operator(void* model, JacobianOperator jacobian_operator) {
   assert(model != nullptr);
   UserModel* user_model = static_cast<UserModel*>(model);
   user_model->jacobian_operator = jacobian_operator;
}

void uno_set_jacobian_transposed_operator(void* model, JacobianTransposedOperator jacobian_transposed_operator) {
   assert(model != nullptr);
   UserModel* user_model = static_cast<UserModel*>(model);
   user_model->jacobian_transposed_operator = jacobian_transposed_operator;
}

void uno_set_lagrangian_hessian(void* model, int32_t number_hessian_nonzeros, char hessian_triangular_part,
      int32_t* hessian_row_indices, int32_t* hessian_column_indices, Hessian lagrangian_hessian,
      double lagrangian_sign_convention) {
   if (number_hessian_nonzeros <= 0) {
      std::cout << "Please specify a positive number of Lagrangian Hessian nonzeros.\n";
      return;
   }
   if (lagrangian_sign_convention != -1. && lagrangian_sign_convention != 1.) {
      std::cout << "Please specify a Lagrangian sign convention in {-1, 1}.\n";
      return;
   }

   assert(model != nullptr);
   UserModel* user_model = static_cast<UserModel*>(model);
   // make sure that the sign convention is consistent with that of the Hessian operator
   if (user_model->lagrangian_hessian_operator != nullptr && user_model->lagrangian_sign_convention != lagrangian_sign_convention) {
      std::cout << "Please specify a Lagrangian sign convention consistent with that of the Hessian operator.\n";
      return;
   }
   user_model->number_hessian_nonzeros = number_hessian_nonzeros;
   user_model->hessian_triangular_part = hessian_triangular_part;
   user_model->hessian_row_indices = hessian_row_indices;
   user_model->hessian_column_indices = hessian_column_indices;
   user_model->lagrangian_hessian = lagrangian_hessian;
   user_model->lagrangian_sign_convention = lagrangian_sign_convention;
}

void uno_set_lagrangian_hessian_operator(void* model, HessianOperator lagrangian_hessian_operator,
      double lagrangian_sign_convention) {
   if (lagrangian_sign_convention != -1. && lagrangian_sign_convention != 1.) {
      std::cout << "Please specify a Lagrangian sign convention in {-1, 1}.\n";
      return;
   }

   assert(model != nullptr);
   UserModel* user_model = static_cast<UserModel*>(model);
   // make sure that the sign convention is consistent with that of the Hessian function
   if (user_model->lagrangian_hessian != nullptr && user_model->lagrangian_sign_convention != lagrangian_sign_convention) {
      std::cout << "Please specify a Lagrangian sign convention consistent with that of the Hessian function.\n";
      return;
   }
   user_model->lagrangian_hessian_operator = lagrangian_hessian_operator;
   user_model->lagrangian_sign_convention = lagrangian_sign_convention;
}

void uno_set_user_data(void* model, void* user_data) {
   assert(model != nullptr);
   UserModel* user_model = static_cast<UserModel*>(model);
   user_model->user_data = user_data;
}

void uno_set_initial_primal_iterate(void* model, double* initial_primal_iterate) {
   assert(model != nullptr);
   UserModel* user_model = static_cast<UserModel*>(model);
   user_model->initial_primal_iterate = initial_primal_iterate;
}

void uno_set_initial_dual_iterate(void* model, double* initial_dual_iterate) {
   assert(model != nullptr);
   UserModel* user_model = static_cast<UserModel*>(model);
   user_model->initial_dual_iterate = initial_dual_iterate;
}

void* uno_create_default_options() {
   uno::Options* options = new uno::Options;
   uno::DefaultOptions::load(*options);
   // determine the default solvers based on the available libraries
   const uno::Options subproblem_solvers = uno::DefaultOptions::determine_subproblem_solvers();
   options->set(subproblem_solvers);
   return options;
}

void uno_set_option(void* options, const char* option_name, const char* option_value) {
   assert(options != nullptr);
   uno::Options* uno_options = static_cast<uno::Options*>(options);
   uno_options->set(option_name, option_value);
}

void* uno_create_solver(const void* options) {
   assert(options != nullptr);
   const uno::Options* uno_options = static_cast<const uno::Options*>(options);
   uno::Logger::set_logger(uno_options->get_string("logger"));
   // TODO number of constraints ??
   return new uno::Uno;
}

void uno_optimize(void* solver, const void* model, const void* options) {
   // check the model
   assert(model != nullptr);
   const UserModel* user_model = static_cast<const UserModel*>(model);
   if (!user_model->objective_function && !user_model->constraint_functions) {
      std::cout << "Please specify at least an objective or constraints.\n";
      return;
   }

   // generate the initial primal-dual iterate
   uno::Iterate initial_iterate(static_cast<size_t>(user_model->number_variables),
      static_cast<size_t>(user_model->number_constraints));
   if (user_model->initial_primal_iterate != nullptr) {
      std::copy_n(user_model->initial_primal_iterate, user_model->number_variables, initial_iterate.primals.begin());
   }
   if (user_model->initial_dual_iterate != nullptr) {
      std::copy_n(user_model->initial_dual_iterate, user_model->number_constraints, initial_iterate.multipliers.constraints.begin());
   }

   assert(solver != nullptr);
   uno::Uno* uno_solver = static_cast<uno::Uno*>(solver);
   assert(options != nullptr);
   const uno::Options* uno_options = static_cast<const uno::Options*>(options);
   UnoModel uno_model(*user_model);
   // TODO return result
   uno::Result result = uno_solver->solve(uno_model, initial_iterate, *uno_options);
}

void uno_destroy_model(void* model) {
   assert(model != nullptr);
   delete static_cast<UserModel*>(model);
}

void uno_destroy_options(void* options) {
   assert(options != nullptr);
   delete static_cast<uno::Options*>(options);
}

void uno_destroy_solver(void* solver) {
   assert(solver != nullptr);
   delete static_cast<uno::Uno*>(solver);
}