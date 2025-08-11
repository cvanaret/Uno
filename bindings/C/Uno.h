// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_C_API
#define UNO_C_API

#ifdef __cplusplus
extern "C" {
#endif

   // objective_function
   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and
   // stores the objective value of "x" in "objective_value".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   int objective_function(int number_variables, const double* x, double* objective_value, void* user_data);

   // constraints
   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the constraint
   // values at "x" in the vector "constraint_values" of size "number_constraints".
   // - returns an integer that is 0 if the evaluations succeeded, and positive otherwise.
   int constraints(int number_variables, int number_constraints, const double* x, double* constraint_values, void* user_data);

   // objective_gradient
   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the dense objective
   // gradient at "x" in the vector "gradient" of size "number_variables".
   // - returns an integer that is 0 if the evaluations succeeded, and positive otherwise.
   int objective_gradient(int number_variables, const double* x, double* gradient, void* user_data);

   // constraint_jacobian_sparsity
   // - takes as inputs two vectors "row_indices" and "column_indices", allocates memory for them, and stores the row and
   // column indices of the constraint Jacobian entries in coordinate (COO) format.
   void constraint_jacobian_sparsity(int* row_indices, int* column_indices, void* user_data);

   // lagrangian_hessian_sparsity
   // - takes as inputs two vectors "row_indices" and "column_indices", allocates memory for them, and stores the row and
   // column indices of the Lagrangian Hessian entries in coordinate (COO) format.
   void lagrangian_hessian_sparsity(int* row_indices, int* column_indices, void* user_data);

   // constraint_jacobian
   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the entries of the
   // sparse constraint Jacobian in the vector "jacobian" of size "number_jacobian_nonzeros". The values should be in
   // the same order as the indices provided in "constraint_jacobian_sparsity".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   int constraint_jacobian(int number_variables, int number_jacobian_nonzeros, const double* x, double* jacobian,
      void* user_data);

   // lagrangian_hessian
   // - takes as inputs a vector "x" of size "number_variables", an objective multiplier, a vector "multipliers" of
   // Lagrange multipliers of size "number_constraints", and an object "user_data", and stores the entries of the
   // sparse Lagrangian Hessian in the vector "hessian" of size "number_hessian_nonzeros". The values should be in
   // the same order as the indices provided in "constraint_jacobian_sparsity".
   // Only the lower triangular part of the symmetric Lagrangian Hessian should be provided.
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   int lagrangian_hessian(int number_variables, int number_constraints, int number_hessian_nonzeros, const double* x,
      double objective_multiplier, const double* multipliers, double* hessian, void* user_data);

   // constraint_jacobian_vector_product
   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Jacobian should be evaluated at "x" (otherwise, the current constraint Jacobian is used), a vector "vector"
   // of size "number_variables" and an object "user_data", and stores the Jacobian-vector product in the vector
   // "result" of size "number_constraints".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   int constraint_jacobian_vector_product(int number_variables, int number_constraints, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* user_data);

   // constraint_jacobian_transposed_vector_product
   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Jacobian should be evaluated at "x" (otherwise, the current constraint Jacobian is used), a vector "vector"
   // of size "number_constraints" and an object "user_data", and stores the Jacobian transposed-vector product in the
   // vector "result" of size "number_variables".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   int constraint_jacobian_transposed_vector_product(int number_variables, int number_constraints, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* user_data);

   // lagrangian_hessian_vector_product
   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Lagrangian Hessian should be evaluated at "x" (otherwise, the current Hessian is used), an objective
   // multiplier, a vector "multipliers" of Lagrange multipliers of size "number_constraints", a vector "vector" of
   // size "number_variables", and an object "user_data", and stores the Hessian-vector product in the vector "result"
   // of size "number_variables".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   int lagrangian_hessian_vector_product(int number_variables, int number_constraints, const double* x,
      bool evaluate_at_x, double objective_multiplier, const double* multipliers, const double* vector,
      double* result, void* user_data);

#ifdef __cplusplus
}
#endif

#endif // UNO_C_API