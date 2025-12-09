@testset "uno -- ADNLPModel" begin
    n = 10
    v = rand(n)
    obj(t) = sum(abs2, t)/2 + dot(t,v)

    nlp = ADNLPModel(obj, ones(n); backend=:generic)
    (model, solver) = uno(nlp; preset="filtersqp")
    primal_solution = Vector{Float64}(undef, nlp.meta.nvar)
    UnoSolver.uno_get_primal_solution(solver, primal_solution)
end

@testset "uno_model -- uno_solver -- uno_optimize" begin
    nlp = CUTEstModel{Float64}("HS15")
    model = uno_model(nlp)
    solver = uno_solver(preset="filtersqp", print_solution=true)
    uno_optimize(solver, model)

    optimization_status = UnoSolver.uno_get_optimization_status(solver)
    solution_status = UnoSolver.uno_get_solution_status(solver)
    solution_objective = UnoSolver.uno_get_solution_objective(solver)
    solution_primal_feasibility = UnoSolver.uno_get_solution_primal_feasibility(solver)
    solution_stationarity = UnoSolver.uno_get_solution_stationarity(solver)
    solution_complementarity = UnoSolver.uno_get_solution_complementarity(solver)

    primal_solution = Vector{Float64}(undef, nlp.meta.nvar)
    UnoSolver.uno_get_primal_solution(solver, primal_solution)

    constraint_dual_solution = Vector{Float64}(undef, nlp.meta.ncon)
    UnoSolver.uno_get_constraint_dual_solution(solver, constraint_dual_solution)

    lower_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
    UnoSolver.uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)

    upper_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
    UnoSolver.uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)

    for i in 1:nlp.meta.nvar
        @test primal_solution[i] == UnoSolver.uno_get_primal_solution_component(solver, i-1)
        @test lower_bound_dual_solution[i] == UnoSolver.uno_get_lower_bound_dual_solution_component(solver, i-1)
        @test upper_bound_dual_solution[i] == UnoSolver.uno_get_upper_bound_dual_solution_component(solver, i-1)
    end

    for j in 1:nlp.meta.ncon
        @test constraint_dual_solution[j] == UnoSolver.uno_get_constraint_dual_solution_component(solver, j-1)
    end
end

@testset "uno -- CUTEst" begin
    nlp = CUTEstModel{Float64}("BYRDSPHR")
    model, solver = uno(nlp, preset="filtersqp", print_solution=true)

    optimization_status = UnoSolver.uno_get_optimization_status(solver)
    solution_status = UnoSolver.uno_get_solution_status(solver)
    solution_objective = UnoSolver.uno_get_solution_objective(solver)
    solution_primal_feasibility = UnoSolver.uno_get_solution_primal_feasibility(solver)
    solution_stationarity = UnoSolver.uno_get_solution_stationarity(solver)
    solution_complementarity = UnoSolver.uno_get_solution_complementarity(solver)

    primal_solution = Vector{Float64}(undef, nlp.meta.nvar)
    UnoSolver.uno_get_primal_solution(solver, primal_solution)

    constraint_dual_solution = Vector{Float64}(undef, nlp.meta.ncon)
    UnoSolver.uno_get_constraint_dual_solution(solver, constraint_dual_solution)

    lower_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
    UnoSolver.uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)

    upper_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
    UnoSolver.uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)

    for i in 1:nlp.meta.nvar
        @test primal_solution[i] == UnoSolver.uno_get_primal_solution_component(solver, i-1)
        @test lower_bound_dual_solution[i] == UnoSolver.uno_get_lower_bound_dual_solution_component(solver, i-1)
        @test upper_bound_dual_solution[i] == UnoSolver.uno_get_upper_bound_dual_solution_component(solver, i-1)
    end

    for j in 1:nlp.meta.ncon
        @test constraint_dual_solution[j] == UnoSolver.uno_get_constraint_dual_solution_component(solver, j-1)
    end
end
