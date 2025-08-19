using Uno, NLPModels, CUTEst

nlp = CUTEstModel{Float64}("HS15")
uno_model = uno(nlp)
uno_set_solver_preset(uno_model, "ipopt")
uno_set_solver_option(uno_model, "unbounded_objective_threshold", "-1e15")
uno_optimize(uno_model)
