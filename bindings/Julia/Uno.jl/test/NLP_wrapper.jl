nlp = CUTEstModel{Float64}("HS15")
model = uno_model(nlp)
Uno.uno_set_solver_preset(model, "funnelsqp")
Uno.uno_set_solver_option(model, "print_solution", "yes")
uno_optimize(model)

optimization_status = Uno.uno_get_optimization_status(model)
solution_status = Uno.uno_get_solution_status(model)
solution_objective = Uno.uno_get_solution_objective(model)
solution_primal_feasibility = Uno.uno_get_solution_primal_feasibility(model)
solution_dual_feasibility = Uno.uno_get_solution_dual_feasibility(model)
solution_complementarity = Uno.uno_get_solution_complementarity(model)

primal_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_primal_solution(model, primal_solution)

constraint_dual_solution = Vector{Float64}(undef, nlp.meta.ncon)
Uno.uno_get_constraint_dual_solution(model, constraint_dual_solution)

lower_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_lower_bound_dual_solution(model, lower_bound_dual_solution)

upper_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_upper_bound_dual_solution(model, upper_bound_dual_solution)
