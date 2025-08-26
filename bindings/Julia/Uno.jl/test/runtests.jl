using Uno, Test
using NLPModels, CUTEst

version = Uno.version()
println("The version of Uno is $version.")

include("C_wrapper.jl")
include("NLP_wrapper.jl")
include("MOI_wrapper.jl")