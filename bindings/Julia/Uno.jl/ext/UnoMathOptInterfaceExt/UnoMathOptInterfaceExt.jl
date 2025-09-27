module UnoMathOptInterfaceExt

import Uno
import MathOptInterface as MOI

function __init__()
    setglobal!(Uno, :Optimizer, Optimizer)
    return
end

include("MOI_utils.jl")
include("MOI_wrapper.jl")

end  # module UnoMathOptInterfaceExt
