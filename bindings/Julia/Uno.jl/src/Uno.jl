module Uno

import Uno_jll
import Uno_jll: libuno
import LinearAlgebra
import OpenBLAS32_jll

function __init__()
    config = LinearAlgebra.BLAS.lbt_get_config()
    if !any(lib -> lib.interface == :lp64, config.loaded_libs)
        LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)
    end
    return
end

export uno

include("libuno.jl")
include("C_wrapper.jl")

global Optimizer

end # module Uno