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
   int objective_function(const double* x, int number_variables, double* objective_value, void* user_data);

   // constraints
   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the constraint
   // values at "x" in the vector "constraint_values" of size "number_constraints".
   // - returns an integer that is 0 if the evaluations succeeded, and positive otherwise.
   int constraints(const double* x, int number_variables, double* constraint_values, int number_constraints, void* user_data);

   // objective_gradient
   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the dense objective
   // gradient at "x" in the vector "gradient" of size "number_variables".
   // - returns an integer that is 0 if the evaluations succeeded, and positive otherwise.
   int objective_gradient(const double* x, int number_variables, double* gradient, void* user_data);

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
   // sparse constraint Jacobian in the vector "jacobian" of size "number_nonzeros". The values should be in the same
   // order as the indices provided in "constraint_jacobian_sparsity".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   int constraint_jacobian(const double* x, int number_variables, double* jacobian, int number_nonzeros,
      void* user_data);

   // lagrangian_hessian
   // - takes as inputs a vector "x" of size "number_variables", an objective multiplier, a vector "multipliers" of
   // Lagrange multipliers of size "number_constraints", and an object "user_data", and stores the entries of the
   // sparse Lagrangian Hessian in the vector "hessian" of size "number_nonzeros". The values should be in the same
   // order as the indices provided in "constraint_jacobian_sparsity".
   // Only the lower triangular part of the symmetric Lagrangian Hessian should be provided.
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   int lagrangian_hessian(const double* x, int number_variables, double objective_multiplier, const double* multipliers,
      int number_constraints, double* hessian, int number_nonzeros, void* user_data);

   // constraint_jacobian_vector_product
   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Jacobian should be evaluated at "x" (otherwise, the current constraint Jacobian is used), a vector "vector"
   // of size "number_variables" and an object "user_data", and stores the Jacobian-vector product in the vector
   // "result" of size "number_constraints".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   int constraint_jacobian_vector_product(const double* x, int number_variables, bool evaluate_at_x,
      const double* vector, double* result, int number_constraints, void* user_data);

   // constraint_jacobian_transposed_vector_product
   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Jacobian should be evaluated at "x" (otherwise, the current constraint Jacobian is used), a vector "vector"
   // of size "number_constraints" and an object "user_data", and stores the Jacobian transposed-vector product in the
   // vector "result" of size "number_variables".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   int constraint_jacobian_transposed_vector_product(const double* x, int number_variables, bool evaluate_at_x,
      const double* vector, int number_constraints, double* result, void* user_data);

   // lagrangian_hessian_vector_product
   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Lagrangian Hessian should be evaluated at "x" (otherwise, the current Hessian is used), an objective
   // multiplier, a vector "multipliers" of Lagrange multipliers of size "number_constraints", a vector "vector" of
   // size "number_variables", and an object "user_data", and stores the Hessian-vector product in the vector "result"
   // of size "number_variables".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   int lagrangian_hessian_vector_product(const double* x, int number_variables, bool evaluate_at_x,
      double objective_multiplier, const double* multipliers, int number_constraints, const double* vector,
      double* result, void* user_data);

#ifdef __cplusplus
}
#endif

#endif // UNO_C_API