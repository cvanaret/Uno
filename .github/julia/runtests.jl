# Copyright (c) 2018-2024: Charlie Vanaret and contributors
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
function Optimizer(options = String["logger=SILENT"])
    return AmplNLWriter.Optimizer(Uno_jll.amplexe, options)
end

Optimizer_LP() = Optimizer(["logger=SILENT", "preset=filterslp", "max_iterations=10000"])

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
        Optimizer;
        exclude = [
            # Remove once https://github.com/cvanaret/Uno/issues/39 is fixed
            "005_010",
            # Okay to exclude forever: AmplNLWriter does not support
            # user-defined functions.
            "006_010",
            # Remove once https://github.com/cvanaret/Uno/issues/38 is fixed
            "007_010",
        ],
        primal_target,
    )
    MINLPTests.test_nlp_expr(
        Optimizer_LP;
        exclude = [
            "001_010",  # Local solution
            "003_014",  # Local solution
            "008_010",  # Local solution
            # Remove once https://github.com/cvanaret/Uno/issues/39 is fixed
            "005_010",
            # Okay to exclude forever: AmplNLWriter does not support
            # user-defined functions.
            "006_010",
            # Remove once https://github.com/cvanaret/Uno/issues/38 is fixed
            "007_010",
        ],
        primal_target,
        objective_tol = 1e-4,
        primal_tol = 1e-4,
    )
    # This function tests convex nonlinear programs. Test failures here should
    # never be allowed, because even local NLP solvers should find the global
    # optimum.
    MINLPTests.test_nlp_cvx_expr(Optimizer; primal_target)
    MINLPTests.test_nlp_cvx_expr(
        Optimizer_LP; 
        primal_target,
        objective_tol = 1e-4,
        primal_tol = 1e-4,
        exclude = ["501_011"],  # Iteration limit
    )
end

# This testset runs the full gamut of MOI.Test.runtests. There are a number of
# tests in here with weird edge cases, so a variety of exclusions are expected.
@testset "MathOptInterface.test" begin
    optimizer = MOI.instantiate(
        Optimizer;
        with_cache_type = Float64,
        with_bridge_type = Float64,
    )
    MOI.Test.runtests(
        optimizer,
        MOI.Test.Config(
            # These are pretty loose tolerances so that all tests pass. There
            # are few tests with weird numerics. If tests fail because of
            # tolerances, it might be okay to make these looser, or you could
            # tighten the tolerances used by Uno.
            atol = 1e-4,
            rtol = 1e-4,
            optimal_status = MOI.LOCALLY_SOLVED,
            infeasible_status = MOI.LOCALLY_INFEASIBLE,
            exclude = Any[
                # It's okay to exclude BasisStatus, since AmplNLWriter doesn't
                # support this information, so too with ObjectiveBound, since
                # Uno is a local NLP solver.
                MOI.VariableBasisStatus,
                MOI.ConstraintBasisStatus,
                MOI.ObjectiveBound,
            ],
        );
        exclude = [
            # ==================================================================
            # The following tests are bugs.
            #
            # We should fix issues in Uno, and then try removing these lines.
            #
            # These tests return OTHER_LIMIT instead of LOCALLY_INFEASIBLE. It
            # might be acceptable to leave this as-is, but it would be better to
            # fix.
            r"^test_conic_NormInfinityCone_INFEASIBLE$",
            r"^test_conic_NormOneCone_INFEASIBLE$",
            r"^test_conic_linear_INFEASIBLE$",
            r"^test_conic_linear_INFEASIBLE_2$",
            r"^test_linear_INFEASIBLE$",
            r"^test_linear_INFEASIBLE_2$",
            r"^test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_",
            # ==================================================================
            # The following tests are okay to exclude forever.
            #
            # Uno does not support integrality, and AmplNLWriter has no way of
            # telling if a particular binary supports integers or not (it
            # defaults to assuming yes).
            "Indicator",
            r"[Ii]nteger",
            "Semicontinuous",
            "Semiinteger",
            "SOS1",
            "SOS2",
            "ZeroOne",
            r"^test_cpsat_",
            r"^test_attribute_SolverVersion$",
            r"^test_nonlinear_invalid$",
            r"^test_basic_VectorNonlinearFunction_",
        ],
    )
end
