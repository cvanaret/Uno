## Uno's Fortran interface: how to use Uno from Fortran

Uno's Fortran interface provides access to the full Uno C API through a lightweight wrapper based on `iso_c_binding`.  
It exposes the same model-building, solver configuration, callback mechanisms, and solution inspection routines as the C interface.

The Fortran interface is provided in the file `uno.f90` and must be included in your source code and linked against the Uno library (`libuno`).

An example is available in [`example_uno.f90`](example_uno.f90).

Start by including the Uno Fortran interface:

```fortran
include 'uno.f90'
```

### Using `uno.f90` as a module (optional)

If you prefer to use the Fortran interface as a module, you can wrap it in a simple module like this:

```fortran
module uno
  include 'uno.f90'
end module
```

Then you can `use uno` in your code instead of `include 'uno.f90'`:

```fortran
use uno
```

**Note:** This approach works with most compilers, but `uno.f90` was originally written as an include file, so some compilers may issue warnings. 

### Building an optimization model

Building an optimization model is incremental and starts with the variables:

```fortran
type(c_ptr) :: model
model = uno_create_model(problem_type, number_variables, &
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
character(len=6) :: hessian_model = "exact"//c_null_char

success = uno_set_solver_integer_option(solver, "max_iterations"//c_null_char, max_iterations)
success = uno_set_solver_double_option(solver, "primal_tolerance"//c_null_char, primal_tolerance)
success = uno_set_solver_bool_option(solver, "print_solution"//c_null_char, print_solution)
success = uno_set_solver_string_option(solver, "hessian_model"//c_null_char, hessian_model)
```

Loading options from a file:

```fortran
logical(c_bool) :: success
character(len=8) :: option_file = "uno.opt"//c_null_char

success = uno_load_solver_option_file(solver, option_file)
```

Getting option values:

```fortran
integer(uno_int) :: max_iterations
real(c_double) :: primal_tolerance
logical(c_bool) :: print_solution

max_iterations = uno_get_solver_integer_option(solver, "max_iterations"//c_null_char)
primal_tolerance   = uno_get_solver_double_option(solver, "primal_tolerance"//c_null_char)
print_solution  = uno_get_solver_bool_option(solver, "print_solution"//c_null_char)
```

Setting a preset:

```fortran
logical(c_bool) :: success
success = uno_set_solver_preset(solver, "filtersqp"//c_null_char)
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

* The Fortran interface closely mirrors the C API: most routines correspond one-to-one to their C counterparts.
* All callbacks must follow the C interoperability rules (`bind(C)` and compatible argument types).
* Strings passed to Uno follow the C convention and must be null-terminated (e.g., `"max_iterations"//c_null_char`).

For a complete list of available routines and constants, see [`uno.f90`](uno.f90).
