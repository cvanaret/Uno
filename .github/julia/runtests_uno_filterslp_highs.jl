# Copyright (c) 2018-2025: Charlie Vanaret and contributors
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

Optimizer_Uno_filterslp() = Optimizer(["logger=SILENT", "preset=filterslp", "LP_solver=HiGHS", "max_iterations=10000", "unbounded_objective_threshold=-1e15", "dual_tolerance=1e-6"])

# This testset runs https://github.com/jump-dev/MINLPTests.jl
@testset "MINLPTests" begin
    primal_target = Dict(
        MINLPTests.FEASIBLE_PROBLEM => MOI.FEASIBLE_POINT,
        # If Uno starts writing a .sol file with an infeasible point, change
        # this to `=> MOI.INFEASIBLE_POINT`
        MINLPTests.INFEASIBLE_PROBLEM => MOI.NO_SOLUTION,
    )
    # This function tests (potentially) non-convex nonlinear programs. The tests
    # are meant to be "easy" in the sense that most NLP solvers can find the
    # same global minimum, but a test failure can sometimes be allowed.
    MINLPTests.test_nlp_expr(
        Optimizer_Uno_filterslp;
        exclude = [
            "001_010",  # Local solution
            "003_014",  # Local solution
            "008_010",  # Local solution
            # no gradient at the initial point
            "005_010",
            # Okay to exclude forever: AmplNLWriter does not support
            # user-defined functions.
            "006_010",
            # Remove once https://github.com/cvanaret/Uno/issues/38 is fixed
            "007_010",
        ],
        primal_target,
        objective_tol = 1e-4,
        primal_tol = 1e-3,
    )
    # This function tests convex nonlinear programs. Test failures here should
    # never be allowed, because even local NLP solvers should find the global
    # optimum.
    MINLPTests.test_nlp_cvx_expr(
        Optimizer_Uno_filterslp; 
        primal_target,
        objective_tol = 1e-4,
        primal_tol = 1e-3,
        exclude = ["501_010", "501_011"],  # Iteration limit
    )
end
