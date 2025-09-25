nlp = CUTEstModel{Float64}("HS15")
uno_model = uno(nlp)
Uno.uno_set_solver_preset(uno_model, "funnelsqp")
Uno.uno_set_solver_option(uno_model, "print_solution", "yes")
Uno.uno_optimize(uno_model)

optimization_status = Uno.uno_get_optimization_status(uno_model)
solution_status = Uno.uno_get_solution_status(uno_model)
solution_objective = Uno.uno_get_solution_objective(uno_model)
solution_primal_feasibility = Uno.uno_get_solution_primal_feasibility(uno_model)
solution_dual_feasibility = Uno.uno_get_solution_dual_feasibility(uno_model)
solution_complementarity = Uno.uno_get_solution_complementarity(uno_model)

primal_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_primal_solution(uno_model, primal_solution)

constraint_dual_solution = Vector{Float64}(undef, nlp.meta.ncon)
Uno.uno_get_constraint_dual_solution(uno_model, constraint_dual_solution)

lower_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_lower_bound_dual_solution(uno_model, lower_bound_dual_solution)

upper_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_upper_bound_dual_solution(uno_model, upper_bound_dual_solution)
