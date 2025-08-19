using Uno, NLPModels, CUTEst

nlp = CUTEstModel{Float64}("HS10")
uno_model = uno(nlp)
