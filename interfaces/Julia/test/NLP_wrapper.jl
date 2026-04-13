@testset "preset = filtersqp" begin
    @testset "uno_model -- uno_solver -- uno_optimize -- HS15 -- $(hessian_model) Hessian" for hessian_model in ("exact", "LBFGS", "LSR1")
        nlp = CUTEstModel{Float64}("HS15")
        model = uno_model(nlp)
        solver = uno_solver(preset="filtersqp", print_solution=true, hessian_model=hessian_model)
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

        finalize(nlp)
    end

    @testset "uno -- CUTEst -- BYRDSPHR -- $(hessian_model) Hessian" for hessian_model in ("exact", "LBFGS", "LSR1")
        nlp = CUTEstModel{Float64}("BYRDSPHR")
        stats = uno(nlp, preset="filtersqp", print_solution=true, hessian_model=hessian_model)
        model, solver = stats.model, stats.solver

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
            @test stats.solution[i] == UnoSolver.uno_get_primal_solution_component(solver, i-1)
            @test stats.multipliers_L[i] == UnoSolver.uno_get_lower_bound_dual_solution_component(solver, i-1)
            @test stats.multipliers_U[i] == UnoSolver.uno_get_upper_bound_dual_solution_component(solver, i-1)
        end

        for j in 1:nlp.meta.ncon
            @test stats.multipliers[j] == UnoSolver.uno_get_constraint_dual_solution_component(solver, j-1)
        end

        finalize(nlp)
    end
end

@testset "preset = ipopt" begin
    for linear_solver in ("MUMPS", "SSIDS", "MA27", "MA57")
        (linear_solver == "SSIDS") && (!haskey(ENV, "OMP_CANCELLATION") || !haskey(ENV, "OMP_PROC_BIND")) && continue
        (linear_solver == "MA27" || linear_solver == "MA57") && !(@ccall HSL_jll.libhsl.LIBHSL_isfunctional()::Bool) && continue

        @testset "linear solver = $linear_solver" begin
            @testset "uno_model -- uno_solver -- uno_optimize -- HS15 -- $(hessian_model) Hessian" for hessian_model in ("exact", "LBFGS")
                nlp = CUTEstModel{Float64}("HS15")
                model = uno_model(nlp)
                solver = uno_solver(preset="ipopt", print_solution=true, linear_solver=linear_solver, hessian_model=hessian_model)
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

                finalize(nlp)
            end

            @testset "uno -- CUTEst -- BYRDSPHR -- $(hessian_model) Hessian" for hessian_model in ("exact", "LBFGS")
                nlp = CUTEstModel{Float64}("BYRDSPHR")
                stats = uno(nlp, preset="ipopt", print_solution=true, linear_solver=linear_solver, hessian_model=hessian_model)
                model, solver = stats.model, stats.solver

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
                    @test stats.solution[i] == UnoSolver.uno_get_primal_solution_component(solver, i-1)
                    @test stats.multipliers_L[i] == UnoSolver.uno_get_lower_bound_dual_solution_component(solver, i-1)
                    @test stats.multipliers_U[i] == UnoSolver.uno_get_upper_bound_dual_solution_component(solver, i-1)
                end

                for j in 1:nlp.meta.ncon
                    @test stats.multipliers[j] == UnoSolver.uno_get_constraint_dual_solution_component(solver, j-1)
                end

                finalize(nlp)
            end
        end
    end
end