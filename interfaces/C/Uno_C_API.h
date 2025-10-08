// Copyright (c) 2025 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_C_API_H
#define UNO_C_API_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
   // Optimization sense
   const int32_t UNO_MINIMIZE =  1;
   const int32_t UNO_MAXIMIZE = -1;

   // Lagrange multiplier sign convention
   const double UNO_MULTIPLIER_POSITIVE =  1.0;
   const double UNO_MULTIPLIER_NEGATIVE = -1.0;

   // Problem type: 'L' = Linear, 'Q' = Quadratic, 'N' = Nonlinear
   const char UNO_PROBLEM_LINEAR    = 'L';
   const char UNO_PROBLEM_QUADRATIC = 'Q';
   const char UNO_PROBLEM_NONLINEAR = 'N';

   // Base indexing style: 0-based (C) or 1-based (Fortran)
   const int32_t UNO_ZERO_BASED_INDEXING = 0;
   const int32_t UNO_ONE_BASED_INDEXING  = 1;

   // Triangular part: 'L' = lower, 'U' = upper
   const char UNO_LOWER_TRIANGLE = 'L';
   const char UNO_UPPER_TRIANGLE = 'U';

   // Optimization status
   const int32_t UNO_SUCCESS = 0;
   const int32_t UNO_ITERATION_LIMIT = 1;
   const int32_t UNO_TIME_LIMIT = 2;
   const int32_t UNO_EVALUATION_ERROR = 3;
   const int32_t UNO_ALGORITHMIC_ERROR = 4;

   // Iterate status
   const int32_t UNO_NOT_OPTIMAL = 0;
   const int32_t UNO_FEASIBLE_KKT_POINT = 1; // feasible stationary point
   const int32_t UNO_FEASIBLE_FJ_POINT = 2; // stationary point without constraint qualification
   const int32_t UNO_INFEASIBLE_STATIONARY_POINT = 3; // infeasible stationary point of constraint violation
   const int32_t UNO_FEASIBLE_SMALL_STEP = 4;
   const int32_t UNO_INFEASIBLE_SMALL_STEP = 5;
   const int32_t UNO_UNBOUNDED = 6;

   // current Uno version is 2.2.1
   const int32_t UNO_VERSION_MAJOR = 2;
   const int32_t UNO_VERSION_MINOR = 2;
   const int32_t UNO_VERSION_PATCH = 1;

   // get the current Uno version as v major.minor.patch
   void uno_get_version(int32_t* major, int32_t* minor, int32_t* patch);

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
   typedef int32_t (*JacobianOperator)(int32_t number_variables, int32_t number_constraints, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Jacobian should be evaluated at "x" (otherwise, the current constraint Jacobian is used), a vector "vector"
   // of size "number_constraints" and an object "user_data", and stores the Jacobian transposed-vector product in the
   // vector "result" of size "number_variables".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int32_t (*JacobianTransposedOperator)(int32_t number_variables, int32_t number_constraints, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Lagrangian Hessian should be evaluated at "x" (otherwise, the current Hessian is used), an objective
   // multiplier, a vector "multipliers" of Lagrange multipliers of size "number_constraints", a vector "vector" of
   // size "number_variables", and an object "user_data", and stores the Hessian-vector product in the vector "result"
   // of size "number_variables".
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise.
   typedef int32_t (*HessianOperator)(int32_t number_variables, int32_t number_constraints, const double* x,
      bool evaluate_at_x, double objective_multiplier, const double* multipliers, const double* vector,
      double* result, void* user_data);

   // creates an optimization model that can be solved by Uno.
   // initially, the model contains "number_variables" variables, no objective function, and no constraints.
   // takes as inputs the type of problem ('L' for linear, 'Q' for quadratic, 'N' for nonlinear), the number of
   // variables, two arrays of lower and upper bounds of size "number_variables", and the vector indexing (0 for
   // C-style indexing, 1 for Fortran-style indexing).
   void* uno_create_model(char problem_type, int32_t number_variables, const double* variables_lower_bounds,
      const double* variables_upper_bounds, int32_t base_indexing);

   // [optional]
   // sets the objective and objective gradient of a given model.
   // takes as inputs the optimization sense (UNO_MINIMIZE or UNO_MAXIMIZE), a function pointer of the objective
   // function and a function pointer of its gradient function.
   // returns true if it succeeded, false otherwise.
   bool uno_set_objective(void* model, int32_t optimization_sense, Objective objective_function,
      ObjectiveGradient objective_gradient);

   // [optional]
   // sets the constraints and constraint Jacobian of a given model.
   // takes as inputs the number of constraints, a function pointer of the constraint functions, two arrays of lower and
   // upper bounds of size "number_constraints", the number of nonzero elements of the Jacobian, two arrays of row and
   // column indices for the constraint Jacobian in COOrdinate format, and a function pointer of the constraint Jacobian.
   // returns true if it succeeded, false otherwise.
   bool uno_set_constraints(void* model, int32_t number_constraints, Constraints constraint_functions,
      const double* constraints_lower_bounds, const double* constraints_upper_bounds, int32_t number_jacobian_nonzeros,
      const int32_t* jacobian_row_indices, const int32_t* jacobian_column_indices, Jacobian constraint_jacobian);

   // [optional]
   // sets the Jacobian operator (computes Jacobian-vector products) of a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_jacobian_operator(void* model, JacobianOperator jacobian_operator);

   // [optional]
   // sets the Jacobian transposed operator (computes Jacobian^T-vector products) of a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_jacobian_transposed_operator(void* model, JacobianTransposedOperator jacobian_transposed_operator);

   // [optional]
   // sets the Lagrangian Hessian of a given model.
   // /!\ since the Lagrangian Hessian is symmetric, we ask for either the lower or the triangular part of the matrix.
   // takes as inputs the number of nonzero elements of the Lagrangian Hessian, a character that specifies whether the
   // lower ('L') or upper ('U') triangular part is provided, two arrays of row and column indices for the Hessian in
   // COOrdinate format, a function pointer of the Hessian, and a scalar in {-1, +1} that determines the sign convention
   // of the Lagrangian:
   // if "lagrangian_sign_convention" == 1,  the Lagrangian is rho*f(x) + y^T c(x)
   // if "lagrangian_sign_convention" == -1, the Lagrangian is rho*f(x) - y^T c(x)
   // returns true if it succeeded, false otherwise.
   bool uno_set_lagrangian_hessian(void* model, int32_t number_hessian_nonzeros, char hessian_triangular_part,
      const int32_t* hessian_row_indices, const int32_t* hessian_column_indices, Hessian lagrangian_hessian,
      double lagrangian_sign_convention);

   // [optional]
   // sets the Lagrangian Hessian operator (computes Hessian-vector products) of a given model.
   // takes as inputs the number of nonzero elements of the Lagrangian Hessian, a function pointer of the Hessian operator,
   // and a scalar in {-1, +1} that determines the sign convention of the Lagrangian:
   // if "lagrangian_sign_convention" == 1,  the Lagrangian is rho*f(x) + y^T c(x)
   // if "lagrangian_sign_convention" == -1, the Lagrangian is rho*f(x) - y^T c(x)
   // returns true if it succeeded, false otherwise.
   bool uno_set_lagrangian_hessian_operator(void* model, int32_t number_hessian_nonzeros,
      HessianOperator lagrangian_hessian_operator, double lagrangian_sign_convention);

   // [optional]
   // sets the user data of a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_user_data(void* model, void* user_data);

   // [optional]
   // sets one component of the initial primal iterate for a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_initial_primal_iterate_component(void* model, int32_t index, double initial_primal_component);

   // [optional]
   // sets one component of the initial dual iterate for a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_initial_dual_iterate_component(void* model, int32_t index, double initial_dual_component);

   // [optional]
   // sets the initial primal iterate of a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_initial_primal_iterate(void* model, const double* initial_primal_iterate);

   // [optional]
   // sets the initial dual iterate of a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_initial_dual_iterate(void* model, const double* initial_dual_iterate);

   // creates the Uno solver.
   void* uno_create_solver();

   // sets a particular option in the Uno solver.
   // takes as inputs the name of the option and the value to which it should be set.
   void uno_set_solver_option(void* solver, const char* option_name, const char* option_value);

   // sets a particular preset in the Uno solver.
   void uno_set_solver_preset(void* solver, const char* preset_name);

   // optimizes a given model using the Uno solver and given options.
   void uno_optimize(void* solver, void* model);

   // gets the optimization status (once the model was solved)
   int32_t uno_get_optimization_status(void* solver);

   // gets the iterate status (once the model was solved)
   int32_t uno_get_solution_status(void* solver);

   // gets the objective value at the solution (once the model was solved)
   double uno_get_solution_objective(void* solver);

   // gets one component of the primal solution (once the model was solved)
   double uno_get_primal_solution_component(void* solver, int32_t index);

   // gets one component of the constraint dual solution (once the model was solved)
   double uno_get_constraint_dual_solution_component(void* solver, int32_t index);

   // gets one component of the lower bound dual solution (once the model was solved)
   double uno_get_lower_bound_dual_solution_component(void* solver, int32_t index);

   // gets one component of the upper bound dual solution (once the model was solved)
   double uno_get_upper_bound_dual_solution_component(void* solver, int32_t index);

   // gets the primal solution (once the model was solved)
   void uno_get_primal_solution(void* solver, double* primal_solution);

   // gets the dual solution associated with the constraints (once the model was solved)
   void uno_get_constraint_dual_solution(void* solver, double* constraint_dual_solution);

   // gets the dual solution associated with the lower bounds (once the model was solved)
   void uno_get_lower_bound_dual_solution(void* solver, double* lower_bound_dual_solution);

   // gets the dual solution associated with the upper bounds (once the model was solved)
   void uno_get_upper_bound_dual_solution(void* solver, double* upper_bound_dual_solution);

   // gets the primal feasibility at the solution (once the model was solved)
   double uno_get_solution_primal_feasibility(void* solver);

   // gets the dual feasibility (aka stationarity) at the solution (once the model was solved)
   double uno_get_solution_dual_feasibility(void* solver);

   // gets the complementarity at the solution (once the model was solved)
   double uno_get_solution_complementarity(void* solver);

   // gets the number of outer iterations required by the solver (once the model was solved)
   int32_t uno_get_number_iterations(void* solver);

   // gets the CPU time required by the solver (once the model was solved)
   double uno_get_cpu_time(void* solver);

   // gets the number of objective evaluations required by the solver (once the model was solved)
   int32_t uno_get_number_objective_evaluations(void* solver);

   // gets the number of constraint evaluations required by the solver (once the model was solved)
   int32_t uno_get_number_constraint_evaluations(void* solver);

   // gets the number of objective gradient evaluations required by the solver (once the model was solved)
   int32_t uno_get_number_objective_gradient_evaluations(void* solver);

   // gets the number of constraint Jacobian evaluations required by the solver (once the model was solved)
   int32_t uno_get_number_jacobian_evaluations(void* solver);

   // gets the number of Lagrangian Hessian evaluations required by the solver (once the model was solved)
   int32_t uno_get_number_hessian_evaluations(void* solver);

   // gets the number of subproblems solved by the solver (once the model was solved)
   int32_t uno_get_number_subproblem_solved_evaluations(void* solver);

   // destroys a given Uno model. Once destroyed, the model cannot be used anymore.
   void uno_destroy_model(void* model);

   // destroy an Uno solver. Once destroyed, the solver cannot be used anymore.
   void uno_destroy_solver(void* solver);

#ifdef __cplusplus
}
#endif

#endif // UNO_C_API_H