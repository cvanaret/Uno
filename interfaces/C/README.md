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

The following optional elements can be added to the model separately:
- the objective function (and its gradient). It is 0 otherwise;
```c
uno_set_objective(model, optimization_sense, objective_function, objective_gradient);
```
- constraint functions (and their Jacobian);
```c
uno_set_constraints(model, number_constraints, constraint_functions,
   constraints_lower_bounds, constraints_upper_bounds, number_jacobian_nonzeros,
   jacobian_row_indices, jacobian_column_indices, constraint_jacobian);
```
- the Lagrangian Hessian;
```c
uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part, 
   hessian_row_indices, hessian_column_indices, lagrangian_hessian, lagrangian_sign_convention);
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
uno_set_lagrangian_hessian_operator(model, number_hessian_nonzeros,
   lagrangian_hessian_operator, lagrangian_sign_convention);
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
uno_set_solver_option(solver, "print_solution", "yes");
```

Setting a preset has Uno mimic an existing solver:
```c
uno_set_solver_preset(solver, "filtersqp");
```

### Solving the model

The model can then be solved by Uno:
```c
uno_optimize(solver, model);
```

### Inspecting the result

A set of functions allows you to inspect the result of the optimization:
- the optimization status (`UNO_SUCCESS`, `UNO_ITERATION_LIMIT`, `UNO_TIME_LIMIT`, `UNO_EVALUATION_ERROR`, `UNO_ALGORITHMIC_ERROR`):
```c
uno_get_optimization_status(solver);
```
- the solution status (`UNO_NOT_OPTIMAL`, `UNO_FEASIBLE_KKT_POINT`, `UNO_FEASIBLE_FJ_POINT`, `UNO_INFEASIBLE_STATIONARY_POINT`, `UNO_FEASIBLE_SMALL_STEP`, `UNO_INFEASIBLE_SMALL_STEP`, `UNO_UNBOUNDED`):
```c
int32_t uno_get_solution_status(solver);
```
- the objective value of the solution:
```c
double uno_get_solution_objective(solver);
```
- the primal solution:
```c
void uno_get_primal_solution(solver, primal_solution);
```
- the dual solution associated with the general constraints:
```c
void uno_get_constraint_dual_solution(solver, constraint_dual_solution);
```
- the dual solution associated with the lower bounds:
```c
void uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution);
```
- the dual solution associated with the upper bounds:
```c
void uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution);
```
- the primal feasibility measure at the solution:
```c
double uno_get_solution_primal_feasibility(solver);
```
- the dual feasibility (aka stationarity) measure at the solution:
```c
double uno_get_solution_dual_feasibility(solver);
```
- the complementarity measure at the solution:
```c
double uno_get_solution_complementarity(solver);
```