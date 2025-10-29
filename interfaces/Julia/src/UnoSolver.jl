module UnoSolver

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

export uno, uno_model, uno_solver, uno_optimize, uno_version
export uno_set_solver_integer_option, uno_set_solver_double_option
export uno_set_solver_bool_option, uno_set_solver_string_option
export uno_set_solver_preset

include("libuno.jl")
include("C_wrapper.jl")

global Optimizer

end  # module UnoSolver
