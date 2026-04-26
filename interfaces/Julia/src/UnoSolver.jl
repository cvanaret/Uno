module UnoSolver

import Uno_jll
import Uno_jll: libuno
import LinearAlgebra
import OpenBLAS32_jll
import SolverCore

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
export uno_set_solver_preset, uno_statistics
export uno_set_initial_primal_iterate, uno_set_initial_dual_iterate
export uno_set_variables_lower_bounds, uno_set_variables_upper_bounds
export uno_set_constraints_lower_bounds, uno_set_constraints_upper_bounds

export UnoExecutionStats

include("libuno.jl")
include("C_wrapper.jl")

global Optimizer

end  # module UnoSolver
