module UnoMathOptInterfaceExt

import Uno
import MathOptInterface as MOI

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

function __init__()
    setglobal!(Uno, :Optimizer, Optimizer)
    return
end

include("MOI_utils.jl")
include("MOI_wrapper.jl")

end  # module UnoMathOptInterfaceExt
