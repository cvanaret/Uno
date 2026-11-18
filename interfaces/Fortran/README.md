## Uno's Fortran interface: how to use Uno from Fortran

Uno provides a Fortran interface built on top of the Uno C API using the module `iso_c_binding`.

To balance completeness and Fortran usability, the interface is split into two files:

* [`uno_c.f90`](uno_c.f90)
* [`uno_fortran.f90`](uno_fortran.f90)

The file `uno_c.f90` defines low-level C bindings for all routines that do not involve C strings, together with all Uno constants, enumerations, and type definitions.
Functions whose signatures contain `char *` arguments or return values are intentionally excluded.

The file `uno_fortran.f90` complements this layer by providing Fortran-friendly wrappers for all string-based C API routines.
These wrappers handle the conversion between Fortran `character` variables and null-terminated C strings, as well as the associated memory management.
Trailing blanks are removed from the strings passed to Uno.
A string that contains `c_null_char` (which C would silently truncate) is rejected by the wrappers without calling the C API: they return the value that the C function returns for an invalid value of that argument (e.g. `UNO_OPTIONS_INVALID_VALUE` for an option value, `UNO_OPTIONS_NOT_FOUND` for an option name, `.false.` for a model name).
A string longer than `UNO_MAX_STRING_LENGTH` (4095) characters is rejected by the C API, with the same status codes.
The strings returned by Uno have at most `UNO_MAX_STRING_LENGTH` characters: the wrappers read at most `UNO_MAX_STRING_LENGTH + 1` characters of them.

Used together, `uno_c.f90` and `uno_fortran.f90` provide complete access to the Uno C API.
Both files must be included when using Uno from Fortran and linked against the Uno library (`libuno`).

An example program is available in [`example_uno.f90`](example_uno.f90).

### Basic usage

Start by including the Fortran interface:

```fortran
include 'uno_c.f90'
...
contains
  include 'uno_fortran.f90'
...
```

### Using the interface as a Fortran module (optional)

You can also wrap the interface in a module for cleaner `use` statements:

```fortran
module uno
  include 'uno_c.f90'
contains
  include 'uno_fortran.f90'
end module
```

Then simply use it in your program:

```fortran
program my_program
  use uno
  ...
end program
```

This approach avoids `include` statements scattered in your code and allows standard module scoping.

**Remark**: We provide the Fortran interface as include files rather than a precompiled module to maximize portability.
Fortran `.mod` files are compiler-specific and can cause issues when cross-compiling or using different compilers.
By including the source directly, users avoid these problems and can build the interface consistently across platforms.

### Building an optimization model

Building an optimization model is incremental and starts with the variables:

```fortran
type(c_ptr) :: model
model = uno_create_model(problem_type, number_variables, &
                         variables_lower_bounds, variables_upper_bounds, &
                         base_indexing)
```
or (for an unconstrained model):
```fortran
model = uno_create_unconstrained_model(problem_type, number_variables, &
                                       base_indexing)
```

The following optional elements can be added or set to the model separately:
* a name for the model (used in Uno's logs and outputs):
```fortran
logical(c_bool) :: success
success = uno_set_model_name(model, name)
```

* lower bounds for the variables:
```fortran
logical(c_bool) :: success
success = uno_set_variables_lower_bounds(model, variables_lower_bounds)
```

* upper bounds for the variables:
```fortran
logical(c_bool) :: success
success = uno_set_variables_upper_bounds(model, variables_upper_bounds)
```

* a lower bound for a given variable:
```fortran
logical(c_bool) :: success
success = uno_set_variable_lower_bound(model, variable_index, lower_bound)
```

* an upper bound for a given variable:
```fortran
logical(c_bool) :: success
success = uno_set_variable_upper_bound(model, variable_index, upper_bound)
```

* the objective function and its gradient:

```fortran
logical(c_bool) :: success
success = uno_set_objective(model, optimization_sense, &
                            objective_function, objective_gradient)
```

* constraint functions and their Jacobian:

```fortran
logical(c_bool) :: success
success = uno_set_constraints(model, number_constraints, constraint_functions, &
                              constraints_lower_bounds, constraints_upper_bounds, &
                              number_jacobian_nonzeros, jacobian_row_indices, &
                              jacobian_column_indices, jacobian)
```

* lower bounds for the constraints:
```fortran
logical(c_bool) :: success
success = uno_set_constraints_lower_bounds(model, constraints_lower_bounds)
```

* upper bounds for the constraints:
```fortran
logical(c_bool) :: success
success = uno_set_constraints_upper_bounds(model, constraints_upper_bounds)
```

* a lower bound for a given constraint:
```fortran
logical(c_bool) :: success
success = uno_set_constraint_lower_bound(model, constraint_index, lower_bound)
```

* an upper bound for a given constraint:
```fortran
logical(c_bool) :: success
success = uno_set_constraint_upper_bound(model, constraint_index, upper_bound)
```

* the Lagrangian Hessian:

```fortran
logical(c_bool) :: success
success = uno_set_lagrangian_hessian(model, number_hessian_nonzeros, &
                                     hessian_triangular_part, &
                                     hessian_row_indices, hessian_column_indices, &
                                     lagrangian_hessian)
```

* a Jacobian operator:

```fortran
logical(c_bool) :: success
success = uno_set_jacobian_operator(model, jacobian_operator)
```

* a Jacobian-transposed operator:

```fortran
logical(c_bool) :: success
success = uno_set_jacobian_transposed_operator(model, jacobian_transposed_operator)
```

* a Hessian operator:

```fortran
logical(c_bool) :: success
success = uno_set_lagrangian_hessian_operator(model, lagrangian_hessian_operator)
```

* the Lagrangian sign convention (default is `UNO_MULTIPLIER_NEGATIVE`):

```fortran
logical(c_bool) :: success
success = uno_set_lagrangian_sign_convention(model, lagrangian_sign_convention)
```

* user-defined data:

```fortran
logical(c_bool) :: success
success = uno_set_user_data(model, user_data)
```

* an initial primal point:

```fortran
logical(c_bool) :: success
success = uno_set_initial_primal_iterate(model, initial_primal_iterate)
```

* an initial dual point:

```fortran
logical(c_bool) :: success
success = uno_set_initial_dual_iterate(model, initial_dual_iterate)
```

The memory for the model is allocated by Uno and must be freed with:

```fortran
call uno_destroy_model(model)
```

### Creating an instance of the Uno solver

Create a solver instance with:

```fortran
type(c_ptr) :: solver
solver = uno_create_solver()
```

The memory for the solver is allocated by Uno must be freed with:

```fortran
call uno_destroy_solver(solver)
```

### Passing options to the Uno solver

Options can be passed individually:

```fortran
integer(uno_int) :: status
integer(uno_int) :: max_iterations = 1000
real(c_double) :: primal_tolerance = 1.0d-6
logical(c_bool) :: print_solution = .true.
character(len=*) :: hessian_model = "exact"

status = uno_set_solver_integer_option(solver, "max_iterations", max_iterations)
status = uno_set_solver_double_option(solver, "primal_tolerance", primal_tolerance)
status = uno_set_solver_bool_option(solver, "print_solution", print_solution)
status = uno_set_solver_string_option(solver, "hessian_model", hessian_model)
```

The functions that set options return a status code: `UNO_OPTIONS_SUCCESS` if the option was set, otherwise `UNO_OPTIONS_INVALID_ARGUMENT`, `UNO_OPTIONS_NOT_FOUND` (unknown option), `UNO_OPTIONS_WRONG_TYPE`, `UNO_OPTIONS_INVALID_VALUE` (e.g. a misspelled value of a string option), `UNO_OPTIONS_INVALID_FILE` or `UNO_OPTIONS_ERROR`. If the status is not `UNO_OPTIONS_SUCCESS`, the options are left unchanged.

Options can also be loaded from a file that contains an `option value` pair per line (the example file `uno_default.opt` in the root directory contains Uno's default options). The options are set in this particular order. Anything that follows a `#` is treated as a comment and is ignored.

```fortran
integer(uno_int) :: status
character(len=*) :: option_file = "uno.opt"

status = uno_load_solver_option_file(solver, option_file)
```

Getting option values:

```fortran
integer(uno_int) :: status
integer(uno_int) :: max_iterations
real(c_double) :: primal_tolerance
logical(c_bool) :: print_solution
character(:), allocatable :: hessian_model

status = uno_get_solver_integer_option(solver, "max_iterations", max_iterations)
status = uno_get_solver_double_option(solver, "primal_tolerance", primal_tolerance)
status = uno_get_solver_bool_option(solver, "print_solution", print_solution)
status = uno_get_solver_string_option(solver, "hessian_model", hessian_model)
```

The getters return `UNO_OPTIONS_SUCCESS` and write the value if the option has a value, and `UNO_OPTIONS_NOT_SET` if the option exists but has no value, e.g. an option set by the `auto` preset before the solve. After a successful solve, and until the options are modified, the getters return the values used by the solve, including the values picked by the `auto` preset and the options overridden by the solver. Otherwise, they return `UNO_OPTIONS_INVALID_ARGUMENT`, `UNO_OPTIONS_NOT_FOUND`, `UNO_OPTIONS_WRONG_TYPE` or `UNO_OPTIONS_ERROR`. The value is written only if the status is `UNO_OPTIONS_SUCCESS`.

Setting a preset:

```fortran
integer(uno_int) :: status
status = uno_set_solver_preset(solver, "filtersqp")
```

### Solving the model

Solve the optimization problem with:

```fortran
integer(uno_int) :: status
status = uno_optimize(solver, model)
```

`uno_optimize` returns `UNO_OPTIMIZE_SUCCESS` if the optimization was carried out (its outcome is given by the optimization status below), otherwise `UNO_OPTIMIZE_INVALID_ARGUMENT` (the solver or the model is not valid), `UNO_OPTIMIZE_INVALID_MODEL` (neither an objective nor constraints), `UNO_OPTIMIZE_INVALID_OPTIONS` (the options cannot be resolved, e.g. the `auto` preset without any subproblem solver available) or `UNO_OPTIMIZE_ERROR`. An error raised during the optimization itself (e.g. a linear solver that is not available in your build) is reported by the optimization status `UNO_ALGORITHMIC_ERROR`.

### Inspecting the result

The following functions allow you to inspect the solution. They write the value into their last argument and return an `integer(uno_int)` status: `UNO_RESULT_SUCCESS` if the value was written, otherwise `UNO_RESULT_INVALID_ARGUMENT`, `UNO_RESULT_NOT_AVAILABLE` (no optimization yet, or `uno_optimize` failed), `UNO_RESULT_INVALID_INDEX` or `UNO_RESULT_ERROR`. The value is written only if the status is `UNO_RESULT_SUCCESS`.

* optimization status:

```fortran
integer(uno_int) :: status, optimization_status
status = uno_get_optimization_status(solver, optimization_status)
```

* solution status:

```fortran
integer(uno_int) :: solution_status
status = uno_get_solution_status(solver, solution_status)
```

* objective value:

```fortran
real(c_double) :: objective
status = uno_get_solution_objective(solver, objective)
```

* primal solution (or one of its components, with `uno_get_primal_solution_component(solver, index, component)`):

```fortran
real(c_double), dimension(number_variables) :: primal_solution
status = uno_get_primal_solution(solver, primal_solution)
```

* dual solutions:

```fortran
real(c_double), dimension(number_variables) :: lower_bound_dual_solution, upper_bound_dual_solution
real(c_double), dimension(number_constraints) :: constraint_dual_solution

status = uno_get_constraint_dual_solution(solver, constraint_dual_solution)
status = uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)
status = uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)
```

* constraint values:

```fortran
real(c_double), dimension(number_constraints) :: constraint_values

status = uno_get_solution_constraints(solver, constraint_values)
```

* optimality measures:

```fortran
real(c_double) :: primal_feasibility, stationarity, complementarity
status = uno_get_solution_primal_feasibility(solver, primal_feasibility)
status = uno_get_solution_stationarity(solver, stationarity)
status = uno_get_solution_complementarity(solver, complementarity)
```

* solver statistics:

```fortran
integer(uno_int) :: number_iterations
real(c_double) :: cpu_time

status = uno_get_number_iterations(solver, number_iterations)
status = uno_get_cpu_time(solver, cpu_time)
```

* description of the method:

```fortran
character(:), allocatable :: method_description
status = uno_get_method_description(solver, method_description)
```
