## Uno's Fortran interface: how to use Uno from Fortran

Uno provides a Fortran interface built on top of the Uno C API using the module `iso_c_binding`.

To balance completeness and Fortran usability, the interface is split into two files:

* [`uno_c.f90`](uno_c.f90)
* [`uno_fortran.f90`](uno_fortran.f90)

The file `uno_c.f90` defines low-level C bindings for all routines that do not involve C strings, together with all Uno constants, enumerations, and type definitions.
Functions whose signatures contain `char *` arguments or return values are intentionally excluded.

The file `uno_fortran.f90` complements this layer by providing Fortran-friendly wrappers for all string-based C API routines.
These wrappers handle the conversion between Fortran `character` variables and null-terminated C strings, as well as the associated memory management.

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
logical(c_bool) :: success
integer(uno_int) :: max_iterations = 1000
real(c_double) :: primal_tolerance = 1.0d-6
logical(c_bool) :: print_solution = .true.
character(len=*) :: hessian_model = "exact"

success = uno_set_solver_integer_option(solver, "max_iterations", max_iterations)
success = uno_set_solver_double_option(solver, "primal_tolerance", primal_tolerance)
success = uno_set_solver_bool_option(solver, "print_solution", print_solution)
success = uno_set_solver_string_option(solver, "hessian_model", hessian_model)
```

Loading options from a file:

```fortran
logical(c_bool) :: success
character(len=*) :: option_file = "uno.opt"

success = uno_load_solver_option_file(solver, option_file)
```

Getting option values:

```fortran
integer(uno_int) :: max_iterations
real(c_double) :: primal_tolerance
logical(c_bool) :: print_solution

max_iterations = uno_get_solver_integer_option(solver, "max_iterations")
primal_tolerance = uno_get_solver_double_option(solver, "primal_tolerance")
print_solution = uno_get_solver_bool_option(solver, "print_solution")
```

Setting a preset:

```fortran
logical(c_bool) :: success
success = uno_set_solver_preset(solver, "filtersqp")
```

### Solving the model

Solve the optimization problem with:

```fortran
call uno_optimize(solver, model)
```

### Inspecting the result

The following routines allow you to inspect the solution:

* optimization status:

```fortran
integer(uno_int) :: optimization_status
optimization_status = uno_get_optimization_status(solver)
```

* solution status:

```fortran
integer(uno_int) :: solution_status
solution_status = uno_get_solution_status(solver)
```

* objective value:

```fortran
real(c_double) :: objective
objective = uno_get_solution_objective(solver)
```

* primal solution:

```fortran
real(c_double), dimension(number_variables) :: primal_solution
call uno_get_primal_solution(solver, primal_solution)
```

* dual solutions:

```fortran
real(c_double), dimension(number_variables) :: lower_bound_dual_solution, upper_bound_dual_solution
real(c_double), dimension(number_constraints) :: constraint_dual_solution

call uno_get_constraint_dual_solution(solver, constraint_dual_solution)
call uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)
call uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)
```

* optimality measures:

```fortran
real(c_double) :: primal_feasibility, stationarity, complementarity
primal_feasibility = uno_get_solution_primal_feasibility(solver)
stationarity = uno_get_solution_stationarity(solver)
complementarity = uno_get_solution_complementarity(solver)
```

* solver statistics:

```fortran
integer(uno_int) :: number_iterations
real(c_double) :: cpu_time

number_iterations  = uno_get_number_iterations(solver)
cpu_time = uno_get_cpu_time(solver)
```
