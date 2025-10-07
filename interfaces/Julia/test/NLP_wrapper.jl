nlp = CUTEstModel{Float64}("HS15")
model = uno_model(nlp)
solver = uno_solver("funnelsqp", print_solution="yes")
uno_optimize(solver, model)

optimization_status = Uno.uno_get_optimization_status(solver)
solution_status = Uno.uno_get_solution_status(solver)
solution_objective = Uno.uno_get_solution_objective(solver)
solution_primal_feasibility = Uno.uno_get_solution_primal_feasibility(solver)
solution_dual_feasibility = Uno.uno_get_solution_dual_feasibility(solver)
solution_complementarity = Uno.uno_get_solution_complementarity(solver)

primal_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_primal_solution(solver, primal_solution)

constraint_dual_solution = Vector{Float64}(undef, nlp.meta.ncon)
Uno.uno_get_constraint_dual_solution(solver, constraint_dual_solution)

lower_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)

upper_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)

for i in 1:nlp.meta.nvar
    @test primal_solution[i] == Uno.uno_get_primal_solution_component(solver, i-1)
    @test lower_bound_dual_solution[i] == Uno.uno_get_lower_bound_dual_solution_component(solver, i-1)
    @test upper_bound_dual_solution[i] == Uno.uno_get_upper_bound_dual_solution_component(solver, i-1)
end

for j in 1:nlp.meta.ncon
    @test constraint_dual_solution[j] == Uno.uno_get_constraint_dual_solution_component(solver, j-1)
end
