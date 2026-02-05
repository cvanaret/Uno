## Uno's Fortran interface: how to use Uno from Fortran

Uno's Fortran interface provides full access to the Uno C API through a lightweight wrapper based on `iso_c_binding`.
It exposes the same functionality as the C interface, including model creation, solver configuration, callback registration, and solution inspection, with some convenient Fortran wrappers for handling strings.

The interface is split into two files: `uno.f90` and `uno_f.f90`.
Include them in your source code and link against the Uno library (`libuno`).

An example program is available in [`example_uno.f90`](example_uno.f90).

### Basic usage

Start by including the Fortran interface:

```fortran
include 'uno.f90'
...
contains
  include 'uno_f.f90'
...
```

### Using the interface as a Fortran module (optional)

You can also wrap the interface in a module for cleaner `use` statements:

```fortran
module uno
  include 'uno.f90'
contains
  include 'uno_f.f90'
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
model = uno_create_model_f(problem_type, number_variables, &
                           variables_lower_bounds, variables_upper_bounds, &
                           base_indexing)
```

The following optional elements can be added to the model separately:

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

* the Lagrangian Hessian:

```fortran
logical(c_bool) :: success
success = uno_set_lagrangian_hessian(model, number_hessian_nonzeros, &
                                     hessian_triangular_part, &
                                     hessian_row_indices, hessian_column_indices, &
                                     lagrangian_hessian, lagrangian_sign_convention)
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
success = uno_set_lagrangian_hessian_operator(model, &
                                              lagrangian_hessian_operator, &
                                              lagrangian_sign_convention)
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

success = uno_set_solver_integer_option_f(solver, "max_iterations", max_iterations)
success = uno_set_solver_double_option_f(solver, "primal_tolerance", primal_tolerance)
success = uno_set_solver_bool_option_f(solver, "print_solution", print_solution)
success = uno_set_solver_string_option_f(solver, "hessian_model", hessian_model)
```

Loading options from a file:

```fortran
logical(c_bool) :: success
character(len=*) :: option_file = "uno.opt"

success = uno_load_solver_option_file_f(solver, option_file)
```

Getting option values:

```fortran
integer(uno_int) :: max_iterations
real(c_double) :: primal_tolerance
logical(c_bool) :: print_solution

max_iterations = uno_get_solver_integer_option_f(solver, "max_iterations")
primal_tolerance = uno_get_solver_double_option_f(solver, "primal_tolerance")
print_solution = uno_get_solver_bool_option_f(solver, "print_solution")
```

Setting a preset:

```fortran
logical(c_bool) :: success
success = uno_set_solver_preset_f(solver, "filtersqp")
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

### Notes

The Fortran interface in `uno.f90` mirrors the Uno C API, including constants, while `uno_f.f90` provides a Fortran-friendly variant for a subset of routines that handle C `char*` inputs or outputs.

The only difference that necessitates `uno_f.f90` is the handling of strings: C `char*` arguments are converted internally and exposed as standard Fortran `character(len=*)` or allocatable strings.
Null-termination is managed automatically.
Routines using this conversion have the suffix `_f` to distinguish them from the direct C bindings.

For a complete list of available routines and constants, see [`uno.f90`](uno.f90) and [`uno_f.f90`](uno_f.f90).
