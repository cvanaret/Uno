// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_C_API_H
#define UNO_C_API_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
   // current Uno version is 2.0.1
   const int32_t uno_version_major = 2;
   const int32_t uno_version_minor = 0;
   const int32_t uno_version_patch = 1;

   // get the current Uno version as v major.minor.patch
   void uno_get_version(int32_t* major, int32_t* minor, int32_t* patch);

   const double uno_Min = 1.;
   const double uno_Max = -1.;

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and
   // stores the objective value of "x" in "objective_value".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int32_t (*Objective)(int32_t number_variables, const double* x, double* objective_value, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the constraint
   // values at "x" in the vector "constraint_values" of size "number_constraints".
   // - returns an integer that is 0 if the evaluations succeeded, and positive otherwise.
   typedef int32_t (*Constraints)(int32_t number_variables, int32_t number_constraints, const double* x, double* constraint_values,
      void* user_data);

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the dense objective
   // gradient at "x" in the vector "gradient" of size "number_variables".
   // - returns an integer that is 0 if the evaluations succeeded, and positive otherwise.
   typedef int32_t (*ObjectiveGradient)(int32_t number_variables, const double* x, double* gradient, void* user_data);

   // - takes as inputs two vectors "row_indices" and "column_indices", allocates memory for them with size
   // "number_jacobian_nonzeros", and stores the row and column indices of the constraint Jacobian entries in
   // coordinate (COO) format.
   // typedef void (*JacobianSparsity)(int32_t number_jacobian_nonzeros, int32_t* row_indices, int32_t* column_indices, void* user_data);

   // - takes as inputs two vectors "row_indices" and "column_indices", allocates memory for them with size
   // "number_hessian_nonzeros", and stores the row and column indices of the Lagrangian Hessian entries in
   // coordinate (COO) format.
   // typedef void (*HessianSparsity)(int32_t number_hessian_nonzeros, int32_t* row_indices, int32_t* column_indices, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the entries of the
   // sparse constraint Jacobian in the vector "jacobian" of size "number_jacobian_nonzeros". The values should be in
   // the same order as the indices provided in "constraint_jacobian_sparsity".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int32_t (*Jacobian)(int32_t number_variables, int32_t number_jacobian_nonzeros, const double* x, double* jacobian,
      void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", an objective multiplier, a vector "multipliers" of
   // Lagrange multipliers of size "number_constraints", and an object "user_data", and stores the entries of the
   // sparse Lagrangian Hessian in the vector "hessian" of size "number_hessian_nonzeros". The values should be in
   // the same order as the indices provided in "constraint_jacobian_sparsity".
   // Only the lower triangular part of the symmetric Lagrangian Hessian should be provided.
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int32_t (*Hessian)(int32_t number_variables, int32_t number_constraints, int32_t number_hessian_nonzeros,
      const double* x, double objective_multiplier, const double* multipliers, double* hessian, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Jacobian should be evaluated at "x" (otherwise, the current constraint Jacobian is used), a vector "vector"
   // of size "number_variables" and an object "user_data", and stores the Jacobian-vector product in the vector
   // "result" of size "number_constraints".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int32_t (*JacobianVectorProduct)(int32_t number_variables, int32_t number_constraints, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Jacobian should be evaluated at "x" (otherwise, the current constraint Jacobian is used), a vector "vector"
   // of size "number_constraints" and an object "user_data", and stores the Jacobian transposed-vector product in the
   // vector "result" of size "number_variables".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int32_t (*JacobianTransposedVectorProduct)(int32_t number_variables, int32_t number_constraints, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Lagrangian Hessian should be evaluated at "x" (otherwise, the current Hessian is used), an objective
   // multiplier, a vector "multipliers" of Lagrange multipliers of size "number_constraints", a vector "vector" of
   // size "number_variables", and an object "user_data", and stores the Hessian-vector product in the vector "result"
   // of size "number_variables".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int32_t (*HessianVectorProduct)(int32_t number_variables, int32_t number_constraints, const double* x,
      bool evaluate_at_x, double objective_multiplier, const double* multipliers, const double* vector,
      double* result, void* user_data);

   // creates an optimization model that can be solved by Uno.
   // initially, the model contains "number_variables" variables, no objective function, and no constraints.
   // takes as inputs the type of problem ('L' for linear, 'Q' for quadratic, 'N' for nonlinear), the number of
   // variables, two arrays of lower and upper bounds of size "number_variables", and the vector indexing (0 for
   // C-style indexing, 1 for Fortran-style indexing).
   void* uno_create_model(char problem_type, int32_t number_variables, double* variables_lower_bounds,
      double* variables_upper_bounds, int32_t vector_indexing);

   // sets the objective and objective gradient of a given model.
   // takes as inputs the objective sense (uno_Min or uno_Max), a function pointer of the objective function and a
   // function pointer of its gradient function.
   void uno_set_objective(void* model, double objective_sense, Objective objective_function, ObjectiveGradient objective_gradient);

   // sets the constraints and constraint Jacobian of a given model.
   // takes as inputs the number of constraints, a function pointer of the constraint functions, two arrays of lower and
   // upper bounds of size "number_constraints", the number of nonzero elements of the Jacobian, two arrays of row and
   // column indices for the constraint Jacobian in COOrdinate format, and a function pointer of the constraint Jacobian.
   void uno_set_constraints(void* model, int32_t number_constraints, Constraints constraint_functions,
      double* constraints_lower_bounds, double* constraints_upper_bounds, int32_t number_jacobian_nonzeros,
      int32_t* jacobian_row_indices, int32_t* jacobian_column_indices, Jacobian constraint_jacobian);

   // sets the Lagrangian Hessian of a given model.
   // /!\ since the Lagrangian Hessian is symmetric, we ask for either the lower or the triangular part of the matrix.
   // takes as inputs the number of nonzero elements of the Lagrangian Hessian, a character that specifies whether the
   // lower ('L') or upper ('U') triangular part is provided, two arrays of row and column indices for the Hessian in
   // COOrdinate format, and a function pointer of the Hessian.
   void uno_set_lagrangian_hessian(void* model, int32_t number_hessian_nonzeros, char hessian_triangular_part,
      int32_t* hessian_row_indices, int32_t* hessian_column_indices, Hessian lagrangian_hessian);

   // sets the user data of a given model.
   void uno_set_user_data(void* model, void* user_data);

   // destroys a given Uno model. Once destroyed, the model cannot be used anymore.
   void uno_destroy_model(void* model);

#ifdef __cplusplus
}
#endif

#endif // UNO_C_API_H