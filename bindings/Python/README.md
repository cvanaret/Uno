## Uno's Python bindings

Uno's Python bindings allows you to solve an optimization model described by callback functions.
The Python module `unopy` is compiled via the command `make unopy` using pybind11, and is (for the moment) installed via `make install`.
An example is available in the file [example_hs015.c](example/example_hs015.py).

### Building an optimization model

Building an optimization model is incremental and starts with the information about the variables:

```python
model = unopy.Model(problem_type, number_variables,
   variables_lower_bounds, variables_upper_bounds, base_indexing)
```

The following optional elements can be added to the model separately:
- the objective function (and its gradient). It is 0 otherwise;
```python
model.set_objective(optimization_sense, objective_function, objective_gradient)
```
- constraint functions (and their Jacobian);
```python
model.set_constraints(number_constraints, constraint_functions,
   constraints_lower_bounds, constraints_upper_bounds, number_jacobian_nonzeros,
   jacobian_row_indices, jacobian_column_indices, constraint_jacobian)
```
- the Lagrangian Hessian;
```python
model.set_lagrangian_hessian(number_hessian_nonzeros, hessian_triangular_part, 
   hessian_row_indices, hessian_column_indices, lagrangian_hessian, lagrangian_sign_convention)
```
- a Jacobian operator (performs Jacobian-vector products);
```python
model.set_jacobian_operator(jacobian_operator)
```
- a Jacobian-transposed operator (performs Jacobian-transposed-vector products);
```python
model.set_jacobian_transposed_operator(jacobian_transposed_operator)
```
- a Hessian operator (performs Hessian-vector products);
```python
model.set_lagrangian_hessian_operator(number_hessian_nonzeros,
   lagrangian_hessian_operator, lagrangian_sign_convention)
```
- user data of an arbitrary type;
```python
model.set_user_data(user_data)
```
- an initial primal point;
```python
model.set_initial_primal_iterate(initial_primal_iterate)
```
- an initial dual point.
```python
model.set_initial_dual_iterate(initial_dual_iterate)
```

*Each of these functions throws an exception upon failure.*

### Creating an instance of the Uno solver

Create an instance of the Uno solver with a simple function call:
```python
uno_solver = unopy.UnoSolver()
```

### Passing options to the Uno solver

Options can be passed to the Uno solver:
```python
uno_solver.set_option("print_solution", "yes")
```

Setting a preset has Uno mimic an existing solver:
```python
uno_solver.set_preset("filtersqp")
```

### Solving the model

The model can then be solved by Uno:
```python
result = uno_solver.optimize(model)
```

### Inspecting the result

To inspect the result of the optimization, read the attributes of the `result` object:
- the optimization status (`UNO_SUCCESS`, `UNO_ITERATION_LIMIT`, `UNO_TIME_LIMIT`, `UNO_EVALUATION_ERROR`, `UNO_ALGORITHMIC_ERROR`): `result.optimization_status`
- the solution status (`UNO_NOT_OPTIMAL`, `UNO_FEASIBLE_KKT_POINT`, `UNO_FEASIBLE_FJ_POINT`, `UNO_INFEASIBLE_STATIONARY_POINT`, `UNO_FEASIBLE_SMALL_STEP`, `UNO_INFEASIBLE_SMALL_STEP`, `UNO_UNBOUNDED`): `result.solution_status`
- the objective value of the solution: `result.solution_objective`
- the primal solution: `result.primal_solution`
- the dual solution associated with the general constraints: `result.constraint_dual_solution`
- the dual solution associated with the lower bounds: `result.lower_bound_dual_solution`
- the dual solution associated with the upper bounds: `result.upper_bound_dual_solution`
- the primal feasibility measure at the solution: `result.solution_primal_feasibility`
- the dual feasibility (aka stationarity) measure at the solution: `result.solution_dual_feasibility`
- the complementarity measure at the solution: `result.solution_complementarity`
- the number of (outer) iterations: `result.number_iterations`
- the CPU time (in seconds): `result.cpu_time`
- the number of objective evaluations: `result.number_objective_evaluations`
- the number of constraint evaluations: `result.number_constraint_evaluations`
- the number of objective gradient evaluations: `result.number_objective_gradient_evaluations`
- the number of Jacobian evaluations: `result.number_jacobian_evaluations`
- the number of Hessian evaluations: `result.number_hessian_evaluations`
- the number of subproblems solved: `result.number_subproblems_solved`
