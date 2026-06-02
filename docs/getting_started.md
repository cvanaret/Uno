# Getting started

## Building a model and the solver

In [Python](../interfaces/python), [Julia](../interfaces/julia), [C](../interfaces/c), and [Fortran](../interfaces/fortran), building an optimization model is incremental: you may attach bounds, an objective, constraints, and an initial primal-dual iterate to a particular model, as well as a Lagrangian Hessian, a Lagrangian Hessian operator, and a Lagrangian convention when the model is defined manually instead of via a modeling framework.
By default, a model is unconstrained, has no objective, no Hessian, and has the Lagrangian convention similar to $\mathcal{L} = f(x) - y^T c(x) - z^T x$.

The code snippets below present the Python and Julia/JuMP implementation of the `hs015` instance.

=== "Python"

    ```py
    import unopy
    Inf = float("inf")
    
    # hs015.mod
    def objective(x):
        return 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2
        
    def constraints(x, constraint_values):
        constraint_values[:] = [x[0]*x[1], x[0] + x[1]**2]
        
    def objective_gradient(x, gradient):
        gradient[:] = [400.*x[0]**3 - 400.*x[0]*x[1] + 2.*x[0] - 2.,
                       200.*(x[1] - x[0]**2)]
    
    def jacobian(x, jacobian_values):
        jacobian_values[:] = [x[1], 1., x[0], 2.*x[1]]
    
    def lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values):
        hessian_values[:] = [objective_multiplier*(1200*x[0]**2 - 400.*x[1] + 2.),
                            -400.*objective_multiplier*x[0] - multipliers[0],
                            200.*objective_multiplier - 2.*multipliers[1]]
    
    def lagrangian_hessian_operator(x, evaluate_at_x, objective_multiplier,
                                    multipliers, vector, result):
        hessian00 = objective_multiplier*(1200*x[0]**2. - 400.*x[1] + 2.)
        hessian10 = -400.*objective_multiplier*x[0] - multipliers[0]
        hessian11 = 200.*objective_multiplier - 2.*multipliers[1]
        result[:] = [hessian00*vector[0] + hessian10*vector[1],
                    hessian10*vector[0] + hessian11*vector[1]]
    
    if __name__ == '__main__':
        number_variables = 2
        number_constraints = 2
    
        # model creation
        model = unopy.Model(unopy.PROBLEM_NONLINEAR, number_variables,
                            unopy.ZERO_BASED_INDEXING)
        model.set_variables_lower_bounds([-Inf, -Inf])
        model.set_variables_upper_bounds([0.5, Inf])
        model.set_objective(unopy.MINIMIZE, objective, objective_gradient)
        model.set_constraints(number_constraints, constraints, [1., 0.], [Inf, Inf],
                              4, [0, 1, 0, 1], [0, 0, 1, 1], jacobian)
        model.set_lagrangian_hessian(3, unopy.LOWER_TRIANGLE, [0, 1, 1], [0, 0, 1],
                                     lagrangian_hessian)
        model.set_lagrangian_sign_convention(unopy.MULTIPLIER_NEGATIVE)
        model.set_initial_primal_iterate([-2., 1.])
    
        # solver creation
        uno_solver = unopy.UnoSolver()
        uno_solver.set_preset("filtersqp")
    
        # solve with the filtersqp preset
        print("Solving with Uno", unopy.current_uno_version())
        result = uno_solver.optimize(model)
        print("Objective at solution:", result.solution_objective)
    ```

=== "Julia/JuMP"

    ```julia
    using UnoSolver, JuMP
    
    jump_model = Model(() -> UnoSolver.Optimizer(preset="filtersqp"))
    x0 = [-2, 1]
    uvar = [0.5, Inf]
    @variable(jump_model, x[i = 1:2] ≤ uvar[i], start = x0[i])
    @objective(jump_model, Min, 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2)
    @constraint(jump_model, x[1] * x[2] - 1 ≥ 0)
    @constraint(jump_model, x[1] + x[2]^2 ≥ 0)
    
    optimize!(jump_model)
    
    termination_status(jump_model)  # solver termination status
    objective_value(jump_model)     # objective value
    value.(x)                       # primal solution
    ```

## Setting options

