# Getting started

The code snippets below present the Python and JuMP implementations of the `hs015` instance.

## Importing the package

=== "Python"

    ```py
    import unopy
    ```

=== "Julia/JuMP"

    ```julia
    using UnoSolver, JuMP
    ```

## Building the solver

The solver is created as follows:

=== "Python"

    ```py
    solver = unopy.UnoSolver()
    ```

=== "Julia/JuMP"

    ```julia
    solver = UnoSolver.Optimizer
    ```

### Setting a preset

Presets are strategy combinations that mimic existing solvers (like filterSQP or IPOPT). A preset may be set as follows:

=== "Python"

    ```py
    solver.set_preset("filtersqp")
    ```

=== "Julia/JuMP"

    ```julia
    solver = () -> UnoSolver.Optimizer(preset="filtersqp")
    ```

### Setting options

Additional options may be set with the following syntax:

=== "Python"

    ```py
    solver.set_options("hessian_model", "LBFGS")
    ```

=== "Julia/JuMP"

    ```julia
    solver = () -> UnoSolver.Optimizer(preset="filtersqp", hessian_model="LBFGS")
    ```

## Building a model

In [Python](../interfaces/python), [Julia](../interfaces/julia), [C](../interfaces/c), and [Fortran](../interfaces/fortran), building an optimization model is incremental: you may attach bounds, an objective, constraints, and an initial primal-dual iterate to a particular model, as well as a Lagrangian Hessian, a Lagrangian Hessian operator, and a Lagrangian convention when the model is defined manually instead of via a modeling framework.
By default, a model is unconstrained, has no objective, no Hessian, and has the Lagrangian convention similar to $\mathcal{L} = f(x) - y^T c(x) - z^T x$.

### Variables and bounds

In Python, the bounds are attached via the functions `set_variables_lower_bounds` and `set_variables_upper_bounds`, while in Julia, the bounds are passed upon the creation of the variables.

=== "Python"

    ```py
    Inf = float("inf")
    number_variables = 2
    model = unopy.Model(unopy.PROBLEM_NONLINEAR, number_variables, unopy.ZERO_BASED_INDEXING)
    model.set_variables_lower_bounds([-Inf, -Inf])
    model.set_variables_upper_bounds([0.5, Inf])
    ```

=== "Julia/JuMP"

    ```julia
    model = Model(solver)
    variables_upper_bounds = [0.5, Inf]
    @variable(model, x[i = 1:2] ≤ variables_upper_bounds[i])
    ```

### Objective

We then attach an objective function. The Python interface expects callback of the objective and its gradient. The derivatives are computed automatically via AD by the JuMP backend.

=== "Python"

    ```py
    def objective(x):
        return 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2

    def objective_gradient(x, gradient):
        gradient[:] = [400.*x[0]**3 - 400.*x[0]*x[1] + 2.*x[0] - 2.,
                       200.*(x[1] - x[0]**2)]
    
    model.set_objective(unopy.MINIMIZE, objective, objective_gradient)
    ```

=== "Julia/JuMP"

    ```julia
    @objective(model, Min, 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2)
    ```

### Constraints

We then attach the constraints and information about the constraint Jacobian, including the sparsity pattern in COO format:

=== "Python"

    ```py
    def constraints(x, constraint_values):
        constraint_values[:] = [x[0]*x[1], x[0] + x[1]**2]

    def jacobian(x, jacobian_values):
        jacobian_values[:] = [x[1], 1., x[0], 2.*x[1]]

    number_constraints = 2
    constraints_lower_bounds = [1., 0.]
    constraints_upper_bounds = [Inf, Inf]
    number_jacobian_nonzeros = 4
    jacobian_row_indices = [0, 1, 0, 1]
    jacobian_column_indices = [0, 0, 1, 1]
    model.set_constraints(number_constraints, constraints, constraints_lower_bounds,
        constraints_upper_bounds, number_jacobian_nonzeros, jacobian_row_indices,
        jacobian_column_indices, jacobian)
    ```

=== "Julia/JuMP"

    ```julia
    @constraint(model, x[1] * x[2] - 1 ≥ 0)
    @constraint(model, x[1] + x[2]^2 ≥ 0)
    ```

### Lagrangian Hessian

If the Hessian of the Lagrangian is available explicitly, it may be passed to Python as a callback along with its sparsity pattern in COO format:

```py
def lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values):
    hessian_values[:] = [objective_multiplier*(1200*x[0]**2 - 400.*x[1] + 2.),
                         -400.*objective_multiplier*x[0] - multipliers[0],
                         200.*objective_multiplier - 2.*multipliers[1]]

number_hessian_nonzeros = 3
hessian_row_indices = [0, 1, 1]
hessian_column_indices = [0, 0, 1]
model.set_lagrangian_hessian(number_hessian_nonzeros, unopy.LOWER_TRIANGLE,
    hessian_row_indices, hessian_column_indices, lagrangian_hessian)
```

The convention of the Lagrangian should be specified in Python whenever the problem has simple bounds (in which cases the corresponding dual variables will be returned with the correct signs) or general constraints (the convention should match the construction of the Lagrangian Hessian if provided).

```py
model.set_lagrangian_sign_convention(unopy.MULTIPLIER_NEGATIVE)
```

### Initial point

=== "Python"

    ```py
    model.set_initial_primal_iterate([-2., 1.])
    ```

=== "Julia/JuMP"

    ```julia
    set_start_value.(x, [-2., 1.])
    ```

## Solving the model

The model is solved by calling the `optimize` function:

=== "Python"

    ```py
    result = solver.optimize(model)
    ```

=== "Julia/JuMP"

    ```julia
    optimize!(model)
    ```

## Inspecting the result

TODO

=== "Python"

    ```py
    print("Objective at solution:", result.solution_objective)
    ```

=== "Julia/JuMP"

    ```julia
    termination_status(model)  # solver termination status
    objective_value(model)     # objective value
    value.(x)                       # primal solution
    ```
