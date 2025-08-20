using Uno, NLPModels, CUTEst

nlp = CUTEstModel{Float64}("HS15")
uno_model = uno(nlp)
Uno.uno_set_solver_preset(uno_model, "ipopt")
Uno.uno_set_solver_option(uno_model, "unbounded_objective_threshold", "-1e15")
Uno.uno_optimize(uno_model)

optimization_status = Uno.uno_get_optimization_status(uno_model)
solution_status = Uno.uno_get_solution_status(uno_model)
solution_objective = Uno.uno_get_solution_objective(uno_model)
solution_primal_feasibility = Uno.uno_get_solution_primal_feasibility(uno_model)
solution_dual_feasibility = Uno.uno_get_solution_dual_feasibility(uno_model)
solution_complementarity = Uno.uno_get_solution_complementarity(uno_model)

# Should be Julia vectors and not C pointers!
# Just check that we can call the C routines in Julia.
Uno.uno_get_primal_solution(uno_model)
Uno.uno_get_constraint_dual_solution(uno_model)
Uno.uno_get_lower_bound_dual_solution(uno_model)
Uno.uno_get_upper_bound_dual_solution(uno_model)
