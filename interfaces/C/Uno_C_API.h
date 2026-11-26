// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_C_API_H
#define UNO_C_API_H

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include "uno_int.h"

// the functions of the C API do not throw exceptions: their declarations are "noexcept" for C++ callers
#ifdef __cplusplus
#define UNO_NOEXCEPT noexcept
#else
#define UNO_NOEXCEPT
#endif

#ifdef __cplusplus
extern "C" {
#endif
   // Optimization sense
   static const uno_int UNO_MINIMIZE =  1;
   static const uno_int UNO_MAXIMIZE = -1;

   // Lagrange multiplier sign convention
   static const uno_int UNO_MULTIPLIER_POSITIVE =  1;
   static const uno_int UNO_MULTIPLIER_NEGATIVE = -1;

   // Problem type: "LP" = Linear Problem, "Q" = Quadratic Problem, "NLP" = Nonlinear Problem
   static const char UNO_PROBLEM_LINEAR[]    = "LP";
   static const char UNO_PROBLEM_QUADRATIC[] = "QP";
   static const char UNO_PROBLEM_NONLINEAR[] = "NLP";

   // Base indexing style: 0-based (C) or 1-based (Fortran)
   static const uno_int UNO_ZERO_BASED_INDEXING = 0;
   static const uno_int UNO_ONE_BASED_INDEXING  = 1;

   // Triangular part: 'L' = lower, 'U' = upper
   static const char UNO_LOWER_TRIANGLE = 'L';
   static const char UNO_UPPER_TRIANGLE = 'U';

   // Option type: 0 = integer, 1 = double, 2 = boolean, 3 = string, -1 = not found
   static const uno_int UNO_OPTION_TYPE_INTEGER = 0;
   static const uno_int UNO_OPTION_TYPE_DOUBLE = 1;
   static const uno_int UNO_OPTION_TYPE_BOOL = 2;
   static const uno_int UNO_OPTION_TYPE_STRING = 3;
   static const uno_int UNO_OPTION_TYPE_NOT_FOUND = -1;

   // Status codes returned by the functions that deal with options
   static const uno_int UNO_OPTIONS_SUCCESS = 0;
   static const uno_int UNO_OPTIONS_INVALID_ARGUMENT = 1; // null solver, option name, option value or file name
   static const uno_int UNO_OPTIONS_NOT_FOUND = 2; // no option has this name
   static const uno_int UNO_OPTIONS_WRONG_TYPE = 3; // the option exists, but has a different type
   static const uno_int UNO_OPTIONS_INVALID_VALUE = 4; // the value is not valid for this option
   static const uno_int UNO_OPTIONS_INVALID_FILE = 5; // the option file cannot be opened or is malformed
   static const uno_int UNO_OPTIONS_ERROR = 6; // any other error
   static const uno_int UNO_OPTIONS_NOT_SET = 7; // the option has no value: no default value, and set neither by the
   // user nor by the preset. With the "auto" preset, the options set by the preset (e.g. globalization_mechanism) have
   // no value until the solve: uno_optimize picks them depending on the model

   // Status codes returned by the functions that get the result of an optimization
   static const uno_int UNO_RESULT_SUCCESS = 0;
   static const uno_int UNO_RESULT_INVALID_ARGUMENT = 1; // null solver or null pointer to the value
   static const uno_int UNO_RESULT_NOT_AVAILABLE = 2; // no result: no optimization yet, or uno_optimize failed.
   // This is not an error: no message is printed, so this status can be used to check whether a result is available
   static const uno_int UNO_RESULT_INVALID_INDEX = 3; // the index of the component is not valid
   static const uno_int UNO_RESULT_ERROR = 4; // any other error

   // Status codes returned by uno_optimize. UNO_OPTIMIZE_SUCCESS means that the optimization was carried out and a
   // result is available: its outcome (success, iteration limit, evaluation error...) is given by
   // uno_get_optimization_status. The other codes mean that the optimization could not be carried out (no result)
   static const uno_int UNO_OPTIMIZE_SUCCESS = 0;
   static const uno_int UNO_OPTIMIZE_INVALID_ARGUMENT = 1; // null solver or model
   static const uno_int UNO_OPTIMIZE_INVALID_MODEL = 2; // the model has neither an objective nor constraints
   static const uno_int UNO_OPTIMIZE_INVALID_OPTIONS = 3; // the options cannot be resolved (e.g. the "auto" preset
   // without any subproblem solver available)
   static const uno_int UNO_OPTIMIZE_ERROR = 4; // any other error (e.g. out of memory)

   // Optimization status
   static const uno_int UNO_SUCCESS = 0;
   static const uno_int UNO_ITERATION_LIMIT = 1;
   static const uno_int UNO_TIME_LIMIT = 2;
   static const uno_int UNO_EVALUATION_ERROR = 3;
   static const uno_int UNO_ALGORITHMIC_ERROR = 4;
   static const uno_int UNO_USER_TERMINATION = 5;
   static const uno_int UNO_USER_INTERRUPTION = 6;

   // Iterate status
   static const uno_int UNO_NOT_OPTIMAL = 0;
   static const uno_int UNO_FEASIBLE_KKT_POINT = 1; // feasible stationary point
   static const uno_int UNO_FEASIBLE_FJ_POINT = 2; // stationary point without constraint qualification
   static const uno_int UNO_INFEASIBLE_STATIONARY_POINT = 3; // infeasible stationary point of constraint violation
   static const uno_int UNO_FEASIBLE_SMALL_STEP = 4;
   static const uno_int UNO_INFEASIBLE_SMALL_STEP = 5;
   static const uno_int UNO_DIVERGING_ITERATE = 6;
   static const uno_int UNO_UNBOUNDED_OBJECTIVE = 7;

   // current Uno version
   static const uno_int UNO_VERSION_MAJOR = 2;
   static const uno_int UNO_VERSION_MINOR = 9;
   static const uno_int UNO_VERSION_PATCH = 3;

   // options iterator
   typedef void* uno_options_iterator;

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and
   // stores the objective value of "x" in "objective_value".
   // - returns 0 if the evaluation succeeded, a positive integer if it failed (Uno may then try another point), or a
   // negative integer to stop the optimization (UNO_USER_TERMINATION).
   typedef uno_int (*uno_objective_callback)(uno_int number_variables, const double* x, double* objective_value, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the constraint
   // values at "x" in the vector "constraint_values" of size "number_constraints".
   // - returns 0 if the evaluations succeeded, a positive integer if they failed (Uno may then try another point), or a
   // negative integer to stop the optimization (UNO_USER_TERMINATION).
   typedef uno_int (*uno_constraints_callback)(uno_int number_variables, uno_int number_constraints, const double* x, double* constraint_values,
      void* user_data);

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the dense objective
   // gradient at "x" in the vector "gradient" of size "number_variables".
   // - returns 0 if the evaluations succeeded, a positive integer if they failed (Uno may then try another point), or a
   // negative integer to stop the optimization (UNO_USER_TERMINATION).
   typedef uno_int (*uno_objective_gradient_callback)(uno_int number_variables, const double* x, double* gradient, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables" and an object "user_data", and stores the entries of the
   // sparse constraint Jacobian in the vector "jacobian" of size "number_jacobian_nonzeros". The values should be in
   // the same order as the indices provided in "jacobian_sparsity".
   // - returns 0 if the evaluation succeeded, a positive integer if it failed (Uno may then try another point), or a
   // negative integer to stop the optimization (UNO_USER_TERMINATION).
   typedef uno_int (*uno_constraints_jacobian_callback)(uno_int number_variables, uno_int number_jacobian_nonzeros, const double* x,
      double* jacobian_values, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", an objective multiplier, a vector "multipliers" of
   // Lagrange multipliers of size "number_constraints", and an object "user_data", and stores the entries of the
   // sparse Lagrangian Hessian in the vector "hessian" of size "number_hessian_nonzeros". The values should be in
   // the same order as the indices provided in "jacobian_sparsity".
   // Only the lower triangular part of the symmetric Lagrangian Hessian should be provided.
   // - returns 0 if the evaluation succeeded, a positive integer if it failed (Uno may then try another point), or a
   // negative integer to stop the optimization (UNO_USER_TERMINATION).
   typedef uno_int (*uno_lagrangian_hessian_callback)(uno_int number_variables, uno_int number_constraints, uno_int number_hessian_nonzeros,
      const double* x, double objective_multiplier, const double* multipliers, double* hessian_values, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Jacobian should be evaluated at "x" (otherwise, the current constraint Jacobian is used), a vector "vector"
   // of size "number_variables" and an object "user_data", and stores the Jacobian-vector product in the vector
   // "result" of size "number_constraints".
   // - returns 0 if the evaluation succeeded, a positive integer if it failed (Uno may then try another point), or a
   // negative integer to stop the optimization (UNO_USER_TERMINATION).
   typedef uno_int (*uno_constraints_jacobian_operator_callback)(uno_int number_variables, uno_int number_constraints, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Jacobian should be evaluated at "x" (otherwise, the current constraint Jacobian is used), a vector "vector"
   // of size "number_constraints" and an object "user_data", and stores the Jacobian transposed-vector product in the
   // vector "result" of size "number_variables".
   // - returns 0 if the evaluation succeeded, a positive integer if it failed (Uno may then try another point), or a
   // negative integer to stop the optimization (UNO_USER_TERMINATION).
   typedef uno_int (*uno_constraints_jacobian_transposed_operator_callback)(uno_int number_variables, uno_int number_constraints, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* user_data);

   // - takes as inputs a vector "x" of size "number_variables", a boolean "evaluate_at_x" that indicates whether
   // the Lagrangian Hessian should be evaluated at "x" (otherwise, the current Hessian is used), an objective
   // multiplier, a vector "multipliers" of Lagrange multipliers of size "number_constraints", a vector "vector" of
   // size "number_variables", and an object "user_data", and stores the Hessian-vector product in the vector "result"
   // of size "number_variables".
   // - returns 0 if the evaluation succeeded, a positive integer if it failed (Uno may then try another point), or a
   // negative integer to stop the optimization (UNO_USER_TERMINATION).
   typedef uno_int (*uno_lagrangian_hessian_operator_callback)(uno_int number_variables, uno_int number_constraints, const double* x,
      bool evaluate_at_x, double objective_multiplier, const double* multipliers, const double* vector,
      double* result, void* user_data);

   // - takes as inputs the number of variables, the number of constraints, a vector "primals" of size "number_variables",
   // the lower and upper bound multipliers of size "number_variables", a vector "constraint_multipliers" of size
   // "number_constraints", an objective multiplier, the primal feasibility residual, the dual feasibility residual, and
   // the complementarity residual.
   // - the iterate is in the space of the model ("number_variables" and "number_constraints" of the model, no slack or
   // elastic variables), and the multipliers follow the Lagrangian sign convention of the model, as in the result.
   typedef void (*uno_notify_acceptable_iterate_callback)(uno_int number_variables, uno_int number_constraints, const double* primals,
      const double* lower_bound_multipliers, const double* upper_bound_multipliers, const double* constraint_multipliers,
      double objective_multiplier, double primal_feasibility_residual, double stationarity_residual,
      double complementarity_residual, void* user_data);

   // - takes as inputs the number of variables, the number of constraints, a vector "primals" of size "number_variables",
   // the lower and upper bound multipliers of size "number_variables", a vector "constraint_multipliers" of size
   // "number_constraints", an objective multiplier, the primal feasibility residual, the dual feasibility residual, and
   // the complementarity residual.
   // - the iterate is in the space of the model, as for uno_notify_acceptable_iterate_callback.
   // - returns 0 to continue the optimization, and a nonzero value to stop it (UNO_USER_TERMINATION).
   typedef uno_int (*uno_termination_callback)(uno_int number_variables, uno_int number_constraints, const double* primals,
      const double* lower_bound_multipliers, const double* upper_bound_multipliers, const double* constraint_multipliers,
      double objective_multiplier, double primal_feasibility_residual, double stationarity_residual,
      double complementarity_residual, void* user_data);

   // - takes as inputs a vector "buffer" of size "length".
   typedef uno_int (*uno_logger_stream_callback)(const char* buffer, uno_int length, void* user_data);

   // get the current Uno version as v major.minor.patch
   void uno_get_version(uno_int* major, uno_int* minor, uno_int* patch) UNO_NOEXCEPT;

   // creates an optimization model that can be solved by Uno.
   // initially, the model contains "number_variables" variables, no objective function, and no constraints.
   // takes as inputs the type of problem (UNO_PROBLEM_LINEAR, UNO_PROBLEM_QUADRATIC, or UNO_PROBLEM_NONLINEAR), the
   // number of variables, two arrays of lower and upper bounds of size "number_variables", and the vector indexing
   // (0 for C-style indexing, 1 for Fortran-style indexing).
   // bounds (of variables and constraints) whose magnitude is at least 1e20 are infinite.
   void* uno_create_model(const char* problem_type, uno_int number_variables, const double* variables_lower_bounds,
      const double* variables_upper_bounds, uno_int base_indexing) UNO_NOEXCEPT;

   // creates an unconstrained optimization model that can be solved by Uno.
   // initially, the model contains "number_variables" variables, no bounds, no objective function, and no constraints.
   // takes as inputs the type of problem (UNO_PROBLEM_LINEAR, UNO_PROBLEM_QUADRATIC, or UNO_PROBLEM_NONLINEAR), the
   // number of variables, and the vector indexing (0 for C-style indexing, 1 for Fortran-style indexing).
   void* uno_create_unconstrained_model(const char* problem_type, uno_int number_variables, uno_int base_indexing) UNO_NOEXCEPT;

   // [optional]
   // sets the name of a given model (used e.g. in logs and KKT system dumps).
   // an empty name falls back to the interface default.
   // returns true if it succeeded, false otherwise.
   bool uno_set_model_name(void* model, const char* name) UNO_NOEXCEPT;

   // [optional]
   // sets the variables lower bounds of a given model.
   // takes as inputs an array of lower bounds of size "number_variables".
   // returns true if it succeeded, false otherwise.
   bool uno_set_variables_lower_bounds(void* model, const double* variables_lower_bounds) UNO_NOEXCEPT;

   // [optional]
   // sets the variables upper bounds of a given model.
   // takes as inputs an array of upper bounds of size "number_variables".
   // returns true if it succeeded, false otherwise.
   bool uno_set_variables_upper_bounds(void* model, const double* variables_upper_bounds) UNO_NOEXCEPT;

   // [optional]
   // sets the lower bound of a given variable.
   // returns true if it succeeded, false otherwise.
   bool uno_set_variable_lower_bound(void* model, uno_int variable_index, double lower_bound) UNO_NOEXCEPT;

   // [optional]
   // sets the upper bound of a given variable.
   // returns true if it succeeded, false otherwise.
   bool uno_set_variable_upper_bound(void* model, uno_int variable_index, double upper_bound) UNO_NOEXCEPT;

   // [optional]
   // sets the objective and objective gradient of a given model.
   // takes as inputs the optimization sense (UNO_MINIMIZE or UNO_MAXIMIZE), a function pointer of the objective
   // function and a function pointer of its gradient function.
   // returns true if it succeeded, false otherwise.
   bool uno_set_objective(void* model, uno_int optimization_sense, uno_objective_callback objective_function,
      uno_objective_gradient_callback objective_gradient) UNO_NOEXCEPT;

   // [optional]
   // sets the constraints and constraint Jacobian of a given model.
   // takes as inputs the number of constraints, a function pointer of the constraint functions, two arrays of lower and
   // upper bounds of size "number_constraints" (NULL means -inf for the lower bounds, +inf for the upper bounds), the
   // number of nonzero elements of the Jacobian, two arrays of row and column indices for the constraint Jacobian in
   // COOrdinate format, and a function pointer of the constraint Jacobian.
   // returns true if it succeeded, false otherwise.
   bool uno_set_constraints(void* model, uno_int number_constraints, uno_constraints_callback constraint_functions,
      const double* constraints_lower_bounds, const double* constraints_upper_bounds, uno_int number_jacobian_nonzeros,
      const uno_int* jacobian_row_indices, const uno_int* jacobian_column_indices, uno_constraints_jacobian_callback jacobian) UNO_NOEXCEPT;

   // [optional]
   // sets the constraints lower bounds of a given model.
   // takes as inputs an array of lower bounds of size "number_constraints".
   // returns true if it succeeded, false otherwise.
   bool uno_set_constraints_lower_bounds(void* model, const double* constraints_lower_bounds) UNO_NOEXCEPT;

   // [optional]
   // sets the constraints upper bounds of a given model.
   // takes as inputs an array of upper bounds of size "number_constraints".
   // returns true if it succeeded, false otherwise.
   bool uno_set_constraints_upper_bounds(void* model, const double* constraints_upper_bounds) UNO_NOEXCEPT;

   // [optional]
   // sets the lower bound of a given constraint.
   // returns true if it succeeded, false otherwise.
   bool uno_set_constraint_lower_bound(void* model, uno_int constraint_index, double lower_bound) UNO_NOEXCEPT;

   // [optional]
   // sets the upper bound of a given constraint.
   // returns true if it succeeded, false otherwise.
   bool uno_set_constraint_upper_bound(void* model, uno_int constraint_index, double upper_bound) UNO_NOEXCEPT;

   // [optional]
   // sets the Jacobian operator (computes Jacobian-vector products) of a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_jacobian_operator(void* model, uno_constraints_jacobian_operator_callback jacobian_operator) UNO_NOEXCEPT;

   // [optional]
   // sets the Jacobian transposed operator (computes Jacobian^T-vector products) of a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_jacobian_transposed_operator(void* model, uno_constraints_jacobian_transposed_operator_callback jacobian_transposed_operator) UNO_NOEXCEPT;

   // [optional]
   // sets the Lagrangian Hessian of a given model.
   // /!\ since the Lagrangian Hessian is symmetric, we ask for either the lower or the triangular part of the matrix.
   // takes as inputs the number of nonzero elements of the Lagrangian Hessian, a character that specifies whether the
   // lower ('L') or upper ('U') triangular part is provided, two arrays of row and column indices for the Hessian in
   // COOrdinate format, and a function pointer of the Hessian.
   // returns true if it succeeded, false otherwise.
   bool uno_set_lagrangian_hessian(void* model, uno_int number_hessian_nonzeros, char hessian_triangular_part,
      const uno_int* hessian_row_indices, const uno_int* hessian_column_indices, uno_lagrangian_hessian_callback lagrangian_hessian) UNO_NOEXCEPT;

   // [optional]
   // sets the Lagrangian Hessian operator (computes Hessian-vector products) of a given model.
   // takes as inputs a function pointer of the Hessian operator.
   // returns true if it succeeded, false otherwise.
   bool uno_set_lagrangian_hessian_operator(void* model, uno_lagrangian_hessian_operator_callback lagrangian_hessian_operator) UNO_NOEXCEPT;

   // [optional]
   // sets the sign convention of the Lagrangian of a given model:
   // if "lagrangian_sign_convention" == UNO_MULTIPLIER_POSITIVE, the Lagrangian is rho*f(x) + y^T c(x)
   // if "lagrangian_sign_convention" == UNO_MULTIPLIER_NEGATIVE, the Lagrangian is rho*f(x) - y^T c(x)
   // if not set, Uno defaults to UNO_MULTIPLIER_NEGATIVE.
   // returns true if it succeeded, false otherwise.
   bool uno_set_lagrangian_sign_convention(void* model, uno_int lagrangian_sign_convention) UNO_NOEXCEPT;

   // [optional]
   // sets the user data of a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_user_data(void* model, void* user_data) UNO_NOEXCEPT;

   // [optional]
   // sets one component of the initial primal iterate for a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_initial_primal_iterate_component(void* model, uno_int index, double initial_primal_component) UNO_NOEXCEPT;

   // [optional]
   // sets one component of the initial dual iterate for a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_initial_dual_iterate_component(void* model, uno_int index, double initial_dual_component) UNO_NOEXCEPT;

   // [optional]
   // sets the initial primal iterate of a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_initial_primal_iterate(void* model, const double* initial_primal_iterate) UNO_NOEXCEPT;

   // [optional]
   // sets the initial dual iterate of a given model.
   // returns true if it succeeded, false otherwise.
   bool uno_set_initial_dual_iterate(void* model, const double* initial_dual_iterate) UNO_NOEXCEPT;

   // creates the Uno solver.
   void* uno_create_solver() UNO_NOEXCEPT;

   // sets a particular option in the Uno solver.
   // takes as inputs the name of the option and the value to which it should be set.
   // the possible types are integer, double, bool and string. The value of a string option with enumerated values
   // (e.g. logger, hessian_model, linear_solver) is checked against the possible values.
   // returns UNO_OPTIONS_SUCCESS if it succeeded, otherwise UNO_OPTIONS_INVALID_ARGUMENT, UNO_OPTIONS_NOT_FOUND,
   // UNO_OPTIONS_WRONG_TYPE, UNO_OPTIONS_INVALID_VALUE or UNO_OPTIONS_ERROR. If it failed, the options are unchanged.
   uno_int uno_set_solver_integer_option(void* solver, const char* option_name, uno_int option_value) UNO_NOEXCEPT;
   uno_int uno_set_solver_double_option(void* solver, const char* option_name, double option_value) UNO_NOEXCEPT;
   uno_int uno_set_solver_bool_option(void* solver, const char* option_name, bool option_value) UNO_NOEXCEPT;
   uno_int uno_set_solver_string_option(void* solver, const char* option_name, const char* option_value) UNO_NOEXCEPT;

   // gets the type of a particular option in the Uno solver.
   // takes as input the name of the option.
   // returns UNO_OPTION_TYPE_INTEGER, UNO_OPTION_TYPE_DOUBLE, UNO_OPTION_TYPE_BOOL or UNO_OPTION_TYPE_STRING, or
   // UNO_OPTION_TYPE_NOT_FOUND if the option does not exist (or the solver or the option name is null).
   uno_int uno_get_solver_option_type(void* solver, const char* option_name) UNO_NOEXCEPT;

   // iterate over the Uno option types and names
   uno_options_iterator uno_options_begin_iterator() UNO_NOEXCEPT;
   uno_options_iterator uno_options_end_iterator() UNO_NOEXCEPT;
   void uno_options_next_iterator(uno_options_iterator* iterator) UNO_NOEXCEPT;
   const char* uno_options_iterator_name(uno_options_iterator iterator) UNO_NOEXCEPT;
   uno_int uno_options_iterator_type(uno_options_iterator iterator) UNO_NOEXCEPT;
   void uno_options_destroy_iterator(uno_options_iterator* iterator) UNO_NOEXCEPT;

   // [optional] loads the options from a given option file.
   // takes as input the name of the option file.
   // returns UNO_OPTIONS_SUCCESS if it succeeded, otherwise UNO_OPTIONS_INVALID_ARGUMENT, UNO_OPTIONS_INVALID_FILE,
   // UNO_OPTIONS_NOT_FOUND, UNO_OPTIONS_INVALID_VALUE or UNO_OPTIONS_ERROR. If it failed, the options are unchanged.
   uno_int uno_load_solver_option_file(void* solver, const char* file_name) UNO_NOEXCEPT;

   // sets a particular preset in the Uno solver.
   // returns UNO_OPTIONS_SUCCESS if it succeeded, otherwise UNO_OPTIONS_INVALID_ARGUMENT, UNO_OPTIONS_INVALID_VALUE
   // or UNO_OPTIONS_ERROR. If it failed, the options are unchanged.
   uno_int uno_set_solver_preset(void* solver, const char* preset_name) UNO_NOEXCEPT;

   // [optional]
   // sets the user callbacks for solver.
   // returns true if it succeeded, false otherwise.
   bool uno_set_solver_callbacks(void* solver, uno_notify_acceptable_iterate_callback notify_acceptable_iterate_callback,
      uno_termination_callback termination_callback, void* user_data) UNO_NOEXCEPT;

   // [optional]
   // sets the logger stream callback.
   // returns true if it succeeded, false otherwise.
   bool uno_set_logger_stream_callback(uno_logger_stream_callback logger_stream_callback, void* user_data) UNO_NOEXCEPT;

   // [optional]
   // resets the logger stream to the standard output
   // returns true if it succeeded, false otherwise.
   bool uno_reset_logger_stream() UNO_NOEXCEPT;

   // optimizes a given model using the Uno solver and given options.
   // returns UNO_OPTIMIZE_SUCCESS if the optimization was carried out (its outcome is given by
   // uno_get_optimization_status), otherwise UNO_OPTIMIZE_INVALID_ARGUMENT, UNO_OPTIMIZE_INVALID_MODEL,
   // UNO_OPTIMIZE_INVALID_OPTIONS or UNO_OPTIMIZE_ERROR. An error raised during the optimization itself (including the
   // construction of the subproblem solvers, e.g. an unavailable linear solver) is reported by the optimization status
   // UNO_ALGORITHMIC_ERROR. The result of a previous optimization is discarded in all cases.
   uno_int uno_optimize(void* solver, void* model) UNO_NOEXCEPT;

   // gets the description of the optimization method used by Uno to solve the model
   // example: TR Fletcher-filter restoration inequality-constrained SQP method with exact Hessian and no inertia correction
   // like the other result getters, returns UNO_RESULT_SUCCESS if it succeeded, otherwise UNO_RESULT_INVALID_ARGUMENT or
   // UNO_RESULT_NOT_AVAILABLE (no optimization yet, or uno_optimize failed). The string is owned by the solver:
   // it remains valid until the next uno_optimize or uno_destroy_solver.
   uno_int uno_get_method_description(void* solver, const char** method_description) UNO_NOEXCEPT;

   // gets the value of a particular option in the Uno solver.
   // takes as inputs the name of the option and a pointer to the value.
   // the possible types are integer, double, bool and string.
   // the value is the one set by the user, otherwise by the preset, otherwise the default value. After a successful
   // uno_optimize, and until the options are modified, the value is the one used by the solve: this includes the values
   // picked by the "auto" preset and, for a string option, the value overridden by the solver (e.g. hessian_model).
   // returns UNO_OPTIONS_SUCCESS if it succeeded, UNO_OPTIONS_NOT_SET if the option has no value, otherwise
   // UNO_OPTIONS_INVALID_ARGUMENT, UNO_OPTIONS_NOT_FOUND, UNO_OPTIONS_WRONG_TYPE or UNO_OPTIONS_ERROR. The value is
   // written only if the status is UNO_OPTIONS_SUCCESS.
   // the string is owned by the solver: it remains valid until the options of the solver are modified, the next
   // uno_optimize, or uno_destroy_solver.
   uno_int uno_get_solver_integer_option(void* solver, const char* option_name, uno_int* option_value) UNO_NOEXCEPT;
   uno_int uno_get_solver_double_option(void* solver, const char* option_name, double* option_value) UNO_NOEXCEPT;
   uno_int uno_get_solver_bool_option(void* solver, const char* option_name, bool* option_value) UNO_NOEXCEPT;
   uno_int uno_get_solver_string_option(void* solver, const char* option_name, const char** option_value) UNO_NOEXCEPT;

   // the following functions get the result of the last optimization (once the model was solved).
   // they write the value into their last argument and return UNO_RESULT_SUCCESS if it succeeded, otherwise
   // UNO_RESULT_INVALID_ARGUMENT, UNO_RESULT_NOT_AVAILABLE (no optimization yet, or uno_optimize failed),
   // UNO_RESULT_INVALID_INDEX (for the components) or UNO_RESULT_ERROR. The value is written only if the status is
   // UNO_RESULT_SUCCESS. The indices of the components follow the base indexing of the model.

   // gets the optimization status (UNO_SUCCESS, UNO_ITERATION_LIMIT, ...)
   uno_int uno_get_optimization_status(void* solver, uno_int* optimization_status) UNO_NOEXCEPT;

   // gets the solution status (UNO_FEASIBLE_KKT_POINT, ...)
   uno_int uno_get_solution_status(void* solver, uno_int* solution_status) UNO_NOEXCEPT;

   // gets the objective value at the solution
   uno_int uno_get_solution_objective(void* solver, double* solution_objective) UNO_NOEXCEPT;

   // gets one component of the primal solution
   uno_int uno_get_primal_solution_component(void* solver, uno_int index, double* primal_component) UNO_NOEXCEPT;

   // gets one component of the constraint dual solution
   uno_int uno_get_constraint_dual_solution_component(void* solver, uno_int index, double* constraint_dual_component) UNO_NOEXCEPT;

   // gets one component of the lower bound dual solution
   uno_int uno_get_lower_bound_dual_solution_component(void* solver, uno_int index, double* lower_bound_dual_component) UNO_NOEXCEPT;

   // gets one component of the upper bound dual solution
   uno_int uno_get_upper_bound_dual_solution_component(void* solver, uno_int index, double* upper_bound_dual_component) UNO_NOEXCEPT;

   // gets the primal solution (array of size "number_variables")
   uno_int uno_get_primal_solution(void* solver, double* primal_solution) UNO_NOEXCEPT;

   // gets the dual solution associated with the constraints (array of size "number_constraints")
   uno_int uno_get_constraint_dual_solution(void* solver, double* constraint_dual_solution) UNO_NOEXCEPT;

   // gets the dual solution associated with the lower bounds (array of size "number_variables")
   uno_int uno_get_lower_bound_dual_solution(void* solver, double* lower_bound_dual_solution) UNO_NOEXCEPT;

   // gets the dual solution associated with the upper bounds (array of size "number_variables")
   uno_int uno_get_upper_bound_dual_solution(void* solver, double* upper_bound_dual_solution) UNO_NOEXCEPT;

   // gets the values of the constraints at the solution (array of size "number_constraints")
   uno_int uno_get_solution_constraints(void* solver, double* constraint_values) UNO_NOEXCEPT;

   // gets the primal feasibility residual at the solution
   uno_int uno_get_solution_primal_feasibility(void* solver, double* primal_feasibility) UNO_NOEXCEPT;

   // gets the stationarity residual at the solution
   uno_int uno_get_solution_stationarity(void* solver, double* stationarity) UNO_NOEXCEPT;

   // gets the complementarity residual at the solution
   uno_int uno_get_solution_complementarity(void* solver, double* complementarity) UNO_NOEXCEPT;

   // gets the number of outer iterations required by the solver
   uno_int uno_get_number_iterations(void* solver, uno_int* number_iterations) UNO_NOEXCEPT;

   // gets the solve time required by the solver
   uno_int uno_get_cpu_time(void* solver, double* cpu_time) UNO_NOEXCEPT;

   // gets the number of objective evaluations required by the solver
   uno_int uno_get_number_objective_evaluations(void* solver, uno_int* number_evaluations) UNO_NOEXCEPT;

   // gets the number of constraint evaluations required by the solver
   uno_int uno_get_number_constraint_evaluations(void* solver, uno_int* number_evaluations) UNO_NOEXCEPT;

   // gets the number of objective gradient evaluations required by the solver
   uno_int uno_get_number_objective_gradient_evaluations(void* solver, uno_int* number_evaluations) UNO_NOEXCEPT;

   // gets the number of constraint Jacobian evaluations required by the solver
   uno_int uno_get_number_jacobian_evaluations(void* solver, uno_int* number_evaluations) UNO_NOEXCEPT;

   // gets the number of Lagrangian Hessian evaluations required by the solver
   uno_int uno_get_number_hessian_evaluations(void* solver, uno_int* number_evaluations) UNO_NOEXCEPT;

   // gets the number of subproblems solved by the solver
   uno_int uno_get_number_subproblem_solved_evaluations(void* solver, uno_int* number_subproblems) UNO_NOEXCEPT;

   // destroys a given Uno model. Once destroyed, the model cannot be used anymore.
   void uno_destroy_model(void* model) UNO_NOEXCEPT;

   // destroy an Uno solver. Once destroyed, the solver cannot be used anymore.
   void uno_destroy_solver(void* solver) UNO_NOEXCEPT;

#ifdef __cplusplus
}
#endif

#endif // UNO_C_API_H
