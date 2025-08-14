// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_C_API_H
#define UNO_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and
   // stores the objective value of "x" in "objective_value".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int (*Objective)(int number_variables, const double* x, double* objective_value, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the constraint
   // values at "x" in the vector "constraint_values" of size "number_constraints".
   // - returns an integer that is 0 if the evaluations succeeded, and positive otherwise.
   typedef int (*Constraints)(int number_variables, int number_constraints, const double* x, double* constraint_values,
      void* user_data);

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the dense objective
   // gradient at "x" in the vector "gradient" of size "number_variables".
   // - returns an integer that is 0 if the evaluations succeeded, and positive otherwise.
   typedef int (*ObjectiveGradient)(int number_variables, const double* x, double* gradient, void* user_data);

   // - takes as inputs two vectors "row_indices" and "column_indices", allocates memory for them with size
   // "number_jacobian_nonzeros", and stores the row and column indices of the constraint Jacobian entries in
   // coordinate (COO) format.
   typedef void (*JacobianSparsity)(int number_jacobian_nonzeros, int* row_indices, int* column_indices, void* user_data);

   // - takes as inputs two vectors "row_indices" and "column_indices", allocates memory for them with size
   // "number_hessian_nonzeros", and stores the row and column indices of the Lagrangian Hessian entries in
   // coordinate (COO) format.
   typedef void (*HessianSparsity)(int number_hessian_nonzeros, int* row_indices, int* column_indices, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the entries of the
   // sparse constraint Jacobian in the vector "jacobian" of size "number_jacobian_nonzeros". The values should be in
   // the same order as the indices provided in "constraint_jacobian_sparsity".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int (*Jacobian)(int number_variables, int number_jacobian_nonzeros, const double* x, double* jacobian,
      void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", an objective multiplier, a vector "multipliers" of
   // Lagrange multipliers of size "number_constraints", and an object "user_data", and stores the entries of the
   // sparse Lagrangian Hessian in the vector "hessian" of size "number_hessian_nonzeros". The values should be in
   // the same order as the indices provided in "constraint_jacobian_sparsity".
   // Only the lower triangular part of the symmetric Lagrangian Hessian should be provided.
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int (*Hessian)(int number_variables, int number_constraints, int number_hessian_nonzeros, const double* x,
      double objective_multiplier, const double* multipliers, double* hessian, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Jacobian should be evaluated at "x" (otherwise, the current constraint Jacobian is used), a vector "vector"
   // of size "number_variables" and an object "user_data", and stores the Jacobian-vector product in the vector
   // "result" of size "number_constraints".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int (*JacobianVectorProduct)(int number_variables, int number_constraints, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Jacobian should be evaluated at "x" (otherwise, the current constraint Jacobian is used), a vector "vector"
   // of size "number_constraints" and an object "user_data", and stores the Jacobian transposed-vector product in the
   // vector "result" of size "number_variables".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int (*JacobianTransposedVectorProduct)(int number_variables, int number_constraints, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Lagrangian Hessian should be evaluated at "x" (otherwise, the current Hessian is used), an objective
   // multiplier, a vector "multipliers" of Lagrange multipliers of size "number_constraints", a vector "vector" of
   // size "number_variables", and an object "user_data", and stores the Hessian-vector product in the vector "result"
   // of size "number_variables".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int (*HessianVectorProduct)(int number_variables, int number_constraints, const double* x,
      bool evaluate_at_x, double objective_multiplier, const double* multipliers, const double* vector,
      double* result, void* user_data);

   // get the current Uno version as v major.minor.patch
   void uno_get_version(int* major, int* minor, int* patch);

   // creates an optimization model that can be solved by Uno.
   void* uno_create_model();

   // destroys a given Uno model. Once destroyed, the model cannot be used anymore.
   void uno_destroy_model(void* uno_model);

#ifdef __cplusplus
}
#endif

#endif // UNO_C_API_H