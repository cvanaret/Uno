using Uno, Test
using NLPModels, MathOptInterface, CUTEst

version = Uno.version()
println("The version of Uno is $version.")

@testset "C interface" begin
  include("C_wrapper.jl")
end

@testset "Interface for NLPModels.jl" begin
  include("NLP_wrapper.jl")
end

#@testset "Interface for MathOptInterface.jl" begin
#  include("MOI_wrapper.jl")
#end
