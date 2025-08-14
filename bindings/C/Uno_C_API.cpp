// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "Uno_C_API.h"

class CModel {
public:
   CModel(char problem_type, int number_variables, double* variables_lower_bounds, double* variables_upper_bounds,
      int vector_indexing):
         problem_type(problem_type),
         vector_indexing(vector_indexing),
         number_variables(number_variables),
         variables_lower_bounds(variables_lower_bounds),
         variables_upper_bounds(variables_upper_bounds) {
   }

   ~CModel() = default;

   const char problem_type; // 'L' for linear, 'Q' for quadratic, 'N' for nonlinear
   const int vector_indexing; // 0 for C-style indexing, 1 for Fortran-style indexing

   // variables
   const int number_variables;
   double* variables_lower_bounds{nullptr};
   double* variables_upper_bounds{nullptr};

   // objective
   Objective objective_function{nullptr};
   ObjectiveGradient objective_gradient{nullptr};

   // constraints
   int number_constraints{0};
   Constraints constraint_functions{nullptr};
   double* constraints_lower_bounds{nullptr};
   double* constraints_upper_bounds{nullptr};
   int number_jacobian_nonzeros{0};
   JacobianSparsity jacobian_sparsity{nullptr};
   Jacobian constraint_jacobian{nullptr};

   // Hessian
   int number_hessian_nonzeros{0};
   HessianSparsity hessian_sparsity{nullptr};
   Hessian lagrangian_hessian{nullptr};

   void* user_data{nullptr};
};

// current version is 2.0.1
void uno_get_version(int* major, int* minor, int* patch) {
   *major = 2;
   *minor = 0;
   *patch = 1;
}

void* uno_create_model(char problem_type, int number_variables, double* variables_lower_bounds,
      double* variables_upper_bounds, int vector_indexing) {
   return new CModel(problem_type, number_variables, variables_lower_bounds, variables_upper_bounds, vector_indexing);
}

void uno_set_objective(void* model, Objective objective_function, ObjectiveGradient objective_gradient) {
   assert(model != nullptr);
   CModel* c_model = static_cast<CModel*>(model);
   c_model->objective_function = objective_function;
   c_model->objective_gradient = objective_gradient;
}

void uno_set_constraints(void* model, int number_constraints, Constraints constraint_functions,
      double* constraints_lower_bounds, double* constraints_upper_bounds, int number_jacobian_nonzeros,
      JacobianSparsity jacobian_sparsity, Jacobian constraint_jacobian) {
   assert(model != nullptr);
   CModel* c_model = static_cast<CModel*>(model);
   c_model->number_constraints = number_constraints;
   c_model->constraint_functions = constraint_functions;
   c_model->constraints_lower_bounds = constraints_lower_bounds;
   c_model->constraints_upper_bounds = constraints_upper_bounds;
   c_model->number_jacobian_nonzeros = number_jacobian_nonzeros;
   c_model->jacobian_sparsity = jacobian_sparsity;
   c_model->constraint_jacobian = constraint_jacobian;
}

void uno_set_lagrangian_hessian(void* model, int number_hessian_nonzeros, HessianSparsity hessian_sparsity,
      Hessian lagrangian_hessian) {
   assert(model != nullptr);
   CModel* c_model = static_cast<CModel*>(model);
   c_model->number_hessian_nonzeros = number_hessian_nonzeros;
   c_model->hessian_sparsity = hessian_sparsity;
   c_model->lagrangian_hessian = lagrangian_hessian;
}

void uno_set_user_data(void* model, void* user_data) {
   assert(model != nullptr);
   CModel* c_model = static_cast<CModel*>(model);
   c_model->user_data = user_data;
}

void uno_destroy_model(void* model) {
   assert(model != nullptr);
   delete static_cast<CModel*>(model);
}