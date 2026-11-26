## Uno's C interface: how to use Uno from C

Uno's C interface allows you to solve an optimization model described by callback functions.
The file `Uno_C_API.cpp` is compiled as part of `libuno` (static or shared) via `make`, and the header `Uno_C_API.h` is installed via `make install`.
An example is available in the file [example_hs015.c](example/example_hs015.c).

Start by including the Uno header:

```c
#include "Uno_C_API.h"
```

### Building an optimization model

Building an optimization model is incremental and starts with the information about the variables:

```c
void* model = uno_create_model(problem_type, number_variables,
   variables_lower_bounds, variables_upper_bounds, base_indexing);
```
or (for an unconstrained model):
```c
void* model = uno_create_unconstrained_model(problem_type, number_variables, base_indexing);
```

The following optional elements can be added or set to the model separately:
- a name for the model (used in Uno's logs and outputs):
```c
uno_set_model_name(model, name)
```
- lower bounds for the variables:
```c
uno_set_variables_lower_bounds(model, variables_lower_bounds)
```
- upper bounds for the variables:
```c
uno_set_variables_upper_bounds(model, variables_upper_bounds)
```
- a lower bound for a given variable:
```c
uno_set_variable_lower_bound(model, variable_index, lower_bound)
```
- an upper bound for a given variable:
```c
uno_set_variable_upper_bound(model, variable_index, upper_bound)
```
- the objective function (and its gradient). It is 0 otherwise;
```c
uno_set_objective(model, optimization_sense, objective_function, objective_gradient);
```
- constraint functions (and their Jacobian);
```c
uno_set_constraints(model, number_constraints, constraint_functions,
   constraints_lower_bounds, constraints_upper_bounds, number_jacobian_nonzeros,
   jacobian_row_indices, jacobian_column_indices, jacobian);
```
- lower bounds for the constraints:
```c
uno_set_constraints_lower_bounds(model, constraints_lower_bounds)
```
- upper bounds for the constraints:
```c
uno_set_constraints_upper_bounds(model, constraints_upper_bounds)
```
- a lower bound for a given constraint:
```c
uno_set_constraint_lower_bound(model, constraint_index, lower_bound)
```
- an upper bound for a given constraint:
```c
uno_set_constraint_upper_bound(model, constraint_index, upper_bound)
```
- the Lagrangian Hessian;
```c
uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part, 
   hessian_row_indices, hessian_column_indices, lagrangian_hessian);
```
- a Jacobian operator (performs Jacobian-vector products);
```c
uno_set_jacobian_operator(model, jacobian_operator);
```
- a Jacobian-transposed operator (performs Jacobian-transposed-vector products);
```c
uno_set_jacobian_transposed_operator(model, jacobian_transposed_operator);
```
- a Hessian operator (performs Hessian-vector products);
```c
uno_set_lagrangian_hessian_operator(model, lagrangian_hessian_operator);
```
- a Lagrangian sign convention (default is `UNO_MULTIPLIER_NEGATIVE`);
```c
uno_set_lagrangian_sign_convention(model, lagrangian_sign_convention);
```
- user data of an arbitrary type (`void*`);
```c
uno_set_user_data(model, user_data);
```
- an initial primal point;
```c
uno_set_initial_primal_iterate(model, initial_primal_iterate);
```
- an initial dual point.
```c
uno_set_initial_dual_iterate(model, initial_dual_iterate);
```

*Each of these functions returns an integer that is 0 upon success and positive upon failure.*

Each evaluation callback (objective, constraints, gradient, Jacobian, Hessian and operators) returns:
- 0 if the evaluation succeeded;
- a positive integer if the evaluation failed: Uno handles the error (for example, by backtracking to a closer point);
- a negative integer to stop the optimization immediately: the optimization status is then `UNO_USER_TERMINATION`, and the
result is the last accepted iterate.

The memory for the model is allocated by the C interface and must be freed by a call to the function:
```c
uno_destroy_model(model);
```

### Creating an instance of the Uno solver

Create an instance of the Uno solver with a simple function call.
```c
void* solver = uno_create_solver();
```

The memory for the solver is allocated by the C interface and must be freed by a call to the function:
```c
uno_destroy_solver(solver);
```

### Passing options to the Uno solver

Options can be passed to the Uno solver:
```c
uno_set_solver_integer_option(solver, "max_iterations", 1000);
uno_set_solver_double_option(solver, "primal_tolerance", 1.0e-6);
uno_set_solver_bool_option(solver, "print_solution", true);
uno_set_solver_string_option(solver, "hessian_model", "exact");
```

Options can also be loaded from a file that contains an `option value` pair per line (the example file `uno_default.opt` in the root directory contains Uno's default options). The options are set in this particular order. Anything that follows a `#` is treated as a comment and is ignored.
```c
uno_load_solver_option_file(solver, "uno.opt");
```

The option setters and `uno_optimize` never throw. The functions that set options (including `uno_load_solver_option_file` and `uno_set_solver_preset`) return a status code:
- `UNO_OPTIONS_SUCCESS`: the option was set;
- `UNO_OPTIONS_INVALID_ARGUMENT`: the solver, the option name, the option value or the file name is null;
- `UNO_OPTIONS_NOT_FOUND`: no option has this name;
- `UNO_OPTIONS_WRONG_TYPE`: the option exists, but has a different type;
- `UNO_OPTIONS_INVALID_VALUE`: the value is not valid for this option. The value of a string option with enumerated values (e.g. `logger`, `hessian_model`, `linear_solver`, `preset`) must be one of its possible values. Whether a solver (e.g. `MA57`) is available in your build is only checked during the optimization (optimization status `UNO_ALGORITHMIC_ERROR`);
- `UNO_OPTIONS_INVALID_FILE`: the option file cannot be opened or is malformed;
- `UNO_OPTIONS_ERROR`: any other error.

If the status is not `UNO_OPTIONS_SUCCESS`, the options are left unchanged and a warning is printed by the logger.
```c
if (uno_set_solver_string_option(solver, "logger", "SOLENT") != UNO_OPTIONS_SUCCESS) {
   // handle the error
}
```

Getting the value of an option from the Uno solver:
```c
uno_int max_iterations;
uno_int status = uno_get_solver_integer_option(solver, "max_iterations", &max_iterations);
double primal_tolerance;
status = uno_get_solver_double_option(solver, "primal_tolerance", &primal_tolerance);
bool print_solution;
status = uno_get_solver_bool_option(solver, "print_solution", &print_solution);
const char* hessian_model;
status = uno_get_solver_string_option(solver, "hessian_model", &hessian_model);
```
The getters return `UNO_OPTIONS_SUCCESS` and write the value if the option has a value. They return `UNO_OPTIONS_NOT_SET` if the option exists but has no value: some options have no default value and are set by the preset (e.g. `globalization_mechanism`, `hessian_model`). With the `auto` preset, Uno picks them in `uno_optimize` depending on the model, so they have no value before the solve. Otherwise, they return `UNO_OPTIONS_INVALID_ARGUMENT`, `UNO_OPTIONS_NOT_FOUND`, `UNO_OPTIONS_WRONG_TYPE` or `UNO_OPTIONS_ERROR`. The value is written only if the status is `UNO_OPTIONS_SUCCESS`.

The value is the one set by the user, otherwise by the preset, otherwise the default value. After a successful solve, and until the options are modified (or `uno_optimize` is called again), the getters return the values used by the solve: these include the values picked by the `auto` preset (e.g. `globalization_mechanism`) and the string options overridden by the solver (e.g., `hessian_model` may contain `LBFGS` if no Hessian is available).
Note: the user-specified options are persistent between consecutive solves (the possibly overridden options are stored separately).
The string returned by `uno_get_solver_string_option` is owned by the solver: it remains valid until the options of the solver are modified, the next call to `uno_optimize`, or `uno_destroy_solver`.

Setting a preset has Uno mimic an existing solver:
```c
uno_set_solver_preset(solver, "filtersqp");
```

Getting the type of an option from the Uno solver:
```c
uno_int uno_get_solver_option_type(solver, "max_iterations");
```

Iterating over the Uno options:
```c
uno_options_iterator uno_options_begin_iterator();
uno_options_iterator uno_options_end_iterator();
void uno_options_next_iterator(&options_iterator);
const char* uno_options_iterator_name(options_iterator);
uno_int uno_options_iterator_type(options_iterator);
void uno_options_destroy_iterator(&options_iterator);
```

The iterator is automatically destroyed when the end is reached in `uno_options_next_iterator` (no call to `uno_options_destroy_iterator` is required). After destruction, the iterator is set to null.

Example of iterating over Uno options:
```c
for (uno_options_iterator it = uno_options_begin_iterator(); it != uno_options_end_iterator(); uno_options_next_iterator(&it)) {
   const char* name = uno_options_iterator_name(it);
   uno_int type = uno_options_iterator_type(it);
   if (type == UNO_OPTION_TYPE_DOUBLE) {
      double value;
      if (uno_get_solver_double_option(solver, name, &value) == UNO_OPTIONS_SUCCESS) {
         // user's code ...
      }
   }
   // ...
}
```

## Setting solver callbacks

Setting the user callbacks to the Uno solver:
```c
uno_set_solver_callbacks(solver, notify_acceptable_iterate_callback, user_termination_callback, user_data);
```
The termination callback is called at each iterate and returns 0 to continue the optimization, or a nonzero value to
stop it (the optimization status is then `UNO_USER_TERMINATION`).
Both callbacks receive the iterate in the space of the model (its numbers of variables and constraints, without the
slack or elastic variables that Uno may introduce), with multipliers in the Lagrangian sign convention of the model, as
returned by the result getters.

Setting the logger stream callback:
```c
uno_set_logger_stream_callback(logger_stream_callback, user_data);
```

and reset the logger stream to the standard output:
```c
uno_reset_logger_stream();
```

### Solving the model

The model can then be solved by Uno:
```c
uno_int status = uno_optimize(solver, model);
```
`uno_optimize` returns a status code:
- `UNO_OPTIMIZE_SUCCESS`: the optimization was carried out and a result is available. Its outcome (success, iteration limit, evaluation error...) is described by the optimization status below;
- `UNO_OPTIMIZE_INVALID_ARGUMENT`: the solver or the model is null;
- `UNO_OPTIMIZE_INVALID_MODEL`: the model has neither an objective nor constraints;
- `UNO_OPTIMIZE_INVALID_OPTIONS`: the options cannot be resolved (e.g. the `auto` preset without any subproblem solver available);
- `UNO_OPTIMIZE_ERROR`: any other error (e.g. out of memory).

In all cases, the result of a previous optimization is discarded. An error raised during the optimization itself, including the construction of the subproblem solvers (e.g. a linear solver that is not available in your build), is reported by the optimization status `UNO_ALGORITHMIC_ERROR`.

### Inspecting the result

A set of functions allows you to inspect the result of the optimization. They write the value into their last argument and return a status code:
- `UNO_RESULT_SUCCESS`: the value was written;
- `UNO_RESULT_INVALID_ARGUMENT`: the solver or the pointer to the value is null;
- `UNO_RESULT_NOT_AVAILABLE`: there is no result (no optimization yet, or `uno_optimize` failed). No message is printed in this case;
- `UNO_RESULT_INVALID_INDEX`: the index of a component is not valid (the indices follow the base indexing of the model);
- `UNO_RESULT_ERROR`: any other error.

The value is written only if the status is `UNO_RESULT_SUCCESS`.

- the optimization status (`UNO_SUCCESS`, `UNO_ITERATION_LIMIT`, `UNO_TIME_LIMIT`, `UNO_EVALUATION_ERROR`, `UNO_ALGORITHMIC_ERROR`, `UNO_USER_TERMINATION`, `UNO_USER_INTERRUPTION`):
```c
uno_int optimization_status;
uno_int status = uno_get_optimization_status(solver, &optimization_status);
```
- the solution status (`UNO_NOT_OPTIMAL`, `UNO_FEASIBLE_KKT_POINT`, `UNO_FEASIBLE_FJ_POINT`, `UNO_INFEASIBLE_STATIONARY_POINT`, `UNO_FEASIBLE_SMALL_STEP`, `UNO_INFEASIBLE_SMALL_STEP`, `UNO_DIVERGING_ITERATE`, `UNO_UNBOUNDED_OBJECTIVE`):
```c
uno_int solution_status;
status = uno_get_solution_status(solver, &solution_status);
```
- the objective value of the solution:
```c
double solution_objective;
status = uno_get_solution_objective(solver, &solution_objective);
```
- the primal solution, or one of its components:
```c
status = uno_get_primal_solution(solver, primal_solution); // array of size number_variables
status = uno_get_primal_solution_component(solver, index, &primal_component);
```
- the dual solution associated with the general constraints, or one of its components:
```c
status = uno_get_constraint_dual_solution(solver, constraint_dual_solution); // array of size number_constraints
status = uno_get_constraint_dual_solution_component(solver, index, &constraint_dual_component);
```
- the dual solution associated with the lower bounds, or one of its components:
```c
status = uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution); // array of size number_variables
status = uno_get_lower_bound_dual_solution_component(solver, index, &lower_bound_dual_component);
```
- the dual solution associated with the upper bounds, or one of its components:
```c
status = uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution); // array of size number_variables
status = uno_get_upper_bound_dual_solution_component(solver, index, &upper_bound_dual_component);
```
- the constraint values at the solution:
```c
status = uno_get_solution_constraints(solver, constraint_values); // array of size number_constraints
```
- the primal feasibility, stationarity and complementarity measures at the solution:
```c
double primal_feasibility, stationarity, complementarity;
status = uno_get_solution_primal_feasibility(solver, &primal_feasibility);
status = uno_get_solution_stationarity(solver, &stationarity);
status = uno_get_solution_complementarity(solver, &complementarity);
```
- the statistics of the solve:
```c
uno_int number_iterations, number_evaluations;
double cpu_time;
status = uno_get_number_iterations(solver, &number_iterations);
status = uno_get_cpu_time(solver, &cpu_time);
status = uno_get_number_objective_evaluations(solver, &number_evaluations);
// also: uno_get_number_constraint_evaluations, uno_get_number_objective_gradient_evaluations,
// uno_get_number_jacobian_evaluations, uno_get_number_hessian_evaluations, uno_get_number_subproblem_solved_evaluations
```
- the description of the method used by Uno (the string is owned by the solver). Like the other getters, it returns `UNO_RESULT_NOT_AVAILABLE` if there is no result:
```c
const char* method_description;
status = uno_get_method_description(solver, &method_description);
```
