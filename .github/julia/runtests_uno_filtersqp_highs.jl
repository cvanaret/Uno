# Copyright (c) 2025: Charlie Vanaret and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# For help with this file, please contact `@odow` on GitHub.

using Test

import AmplNLWriter
import MathOptInterface as MOI
import MINLPTests
import Uno_jll

"""
    Optimizer()

Create a new `AmplNLWriter.Optimizer` object that uses Uno as the backing
solver.
"""
function Optimizer(options)
    return AmplNLWriter.Optimizer(Uno_jll.amplexe, options)
end

# Note: the HiGHS QP solver tolerates only positive semi-definite Hessians.
# We therefore run only on convex instances (for now)
Optimizer_Uno_filtersqp() = Optimizer(["logger=SILENT", "preset=filtersqp", "QP_solver=HiGHS", "max_iterations=10000"])

# This testset runs https://github.com/jump-dev/MINLPTests.jl
@testset "MINLPTests" begin
    primal_target = Dict(
        MINLPTests.FEASIBLE_PROBLEM => MOI.FEASIBLE_POINT,
        MINLPTests.INFEASIBLE_PROBLEM => MOI.INFEASIBLE_POINT,
    )
    primal_tol = 1e-4
    objective_tol = 1e-4
    # This function tests convex nonlinear programs. Test failures here should
    # never be allowed, because even local NLP solvers should find the global
    # optimum.
    MINLPTests.test_nlp_cvx_expr(Optimizer_Uno_filtersqp;
        exclude = [
            "102_010", "105_010", "201_011", "501_010", # negative curvature
            "110_010" # HiGHS cannot solve the first QP (https://github.com/ERGO-Code/HiGHS/issues/2489)
        ],
        primal_target,
        primal_tol,
        objective_tol)
end
