module UnoSolverMathOptInterfaceExt

import UnoSolver
import MathOptInterface as MOI

function __init__()
    setglobal!(UnoSolver, :Optimizer, Optimizer)
    return
end

include("MOI_utils.jl")
include("MOI_wrapper.jl")

end  # module UnoSolverMathOptInterfaceExt
