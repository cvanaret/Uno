// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <iostream>
#include "Uno_C_API.h"
#include "Uno.hpp"
#include "options/DefaultOptions.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"

class CModel {
public:
   CModel(char problem_type, int32_t number_variables, double* variables_lower_bounds, double* variables_upper_bounds,
      int32_t vector_indexing):
         problem_type(problem_type),
         vector_indexing(vector_indexing),
         number_variables(number_variables),
         variables_lower_bounds(variables_lower_bounds),
         variables_upper_bounds(variables_upper_bounds) {
   }

   ~CModel() = default;

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

   // Hessian
   int32_t number_hessian_nonzeros{0};
   // lower ('L') or upper ('U')
   char hessian_triangular_part{}; // default is empty
   int32_t* hessian_row_indices{nullptr};
   int32_t* hessian_column_indices{nullptr};
   Hessian lagrangian_hessian{nullptr};
   double lagrangian_sign_convention{-1.};

   void* user_data{nullptr};
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
   return new CModel(problem_type, number_variables, variables_lower_bounds, variables_upper_bounds, vector_indexing);
}

void uno_set_objective(void* model, double objective_sense, Objective objective_function, ObjectiveGradient objective_gradient) {
   if (objective_sense != 1. && objective_sense != -1.) {
      std::cout << "Please specify a valid objective sense.\n";
      return;
   }

   assert(model != nullptr);
   CModel* c_model = static_cast<CModel*>(model);
   c_model->objective_sense = objective_sense;
   c_model->objective_function = objective_function;
   c_model->objective_gradient = objective_gradient;
}

void uno_set_constraints(void* model, int32_t number_constraints, Constraints constraint_functions,
      double* constraints_lower_bounds, double* constraints_upper_bounds, int32_t number_jacobian_nonzeros,
      int32_t* jacobian_row_indices, int32_t* jacobian_column_indices, Jacobian constraint_jacobian) {
   if (number_constraints <= 0) {
      std::cout << "Please specify a positive number of constraints.\n";
      return;
   }

   assert(model != nullptr);
   CModel* c_model = static_cast<CModel*>(model);
   c_model->number_constraints = number_constraints;
   c_model->constraint_functions = constraint_functions;
   c_model->constraints_lower_bounds = constraints_lower_bounds;
   c_model->constraints_upper_bounds = constraints_upper_bounds;
   c_model->number_jacobian_nonzeros = number_jacobian_nonzeros;
   c_model->jacobian_row_indices = jacobian_row_indices;
   c_model->jacobian_column_indices = jacobian_column_indices;
   c_model->constraint_jacobian = constraint_jacobian;
}

void uno_set_lagrangian_hessian(void* model, int32_t number_hessian_nonzeros, char hessian_triangular_part,
      int32_t* hessian_row_indices, int32_t* hessian_column_indices, Hessian lagrangian_hessian,
      double lagrangian_sign_convention) {
   if (number_hessian_nonzeros <= 0) {
      std::cout << "Please specify a positive number of Lagrangian Hessian nonzeros.\n";
   }
   if (lagrangian_sign_convention != -1. && lagrangian_sign_convention != 1.) {
      std::cout << "Please specify a Lagrangian sign convention in {-1, 1}.\n";
   }

   assert(model != nullptr);
   CModel* c_model = static_cast<CModel*>(model);
   c_model->number_hessian_nonzeros = number_hessian_nonzeros;
   c_model->hessian_triangular_part = hessian_triangular_part;
   c_model->hessian_row_indices = hessian_row_indices;
   c_model->hessian_column_indices = hessian_column_indices;
   c_model->lagrangian_hessian = lagrangian_hessian;
   c_model->lagrangian_sign_convention = lagrangian_sign_convention;
}

void uno_set_user_data(void* model, void* user_data) {
   assert(model != nullptr);
   CModel* c_model = static_cast<CModel*>(model);
   c_model->user_data = user_data;
}

void* uno_create_options() {
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

void uno_optimize(void* solver, const void* model, const void* options, const double* primal_iterate,
      const double* dual_iterate) {
   // check the model
   assert(model != nullptr);
   const CModel* c_model = static_cast<const CModel*>(model);
   if (!c_model->objective_function && !c_model->constraint_functions) {
      std::cout << "Please specify at least an objective or constraints.\n";
      return;
   }

   // generate the initial primal-dual iterate
   uno::Iterate initial_iterate(static_cast<size_t>(c_model->number_variables),
      static_cast<size_t>(c_model->number_constraints));
   if (primal_iterate != nullptr) {
      std::copy_n(primal_iterate, c_model->number_variables, initial_iterate.primals.begin());
   }
   if (dual_iterate != nullptr) {
      std::copy_n(dual_iterate, c_model->number_constraints, initial_iterate.multipliers.constraints.begin());
   }

   assert(solver != nullptr);
   uno::Uno* uno_solver = static_cast<uno::Uno*>(solver);
   assert(options != nullptr);
   const uno::Options* uno_options = static_cast<const uno::Options*>(options);

   // TODO return result
   //uno::Result result = uno_solver->solve(*c_model, initial_iterate, *uno_options);
}

void uno_destroy_model(void* model) {
   assert(model != nullptr);
   delete static_cast<CModel*>(model);
}

void uno_destroy_options(void* options) {
   assert(options != nullptr);
   delete static_cast<uno::Options*>(options);
}

void uno_destroy_solver(void* solver) {
   assert(solver != nullptr);
   delete static_cast<uno::Uno*>(solver);
}