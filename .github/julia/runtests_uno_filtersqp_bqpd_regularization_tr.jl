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

Optimizer_Uno_filtersqp() = Optimizer(["logger=SILENT", "preset=filtersqp", "QP_solver=BQPD", "inertia_correction_strategy=primal", "globalization_mechanism=TR", "max_iterations=10000", "unbounded_objective_threshold=-1e15"])

# This testset runs https://github.com/jump-dev/MINLPTests.jl
@testset "MINLPTests" begin
    primal_target = Dict(
        MINLPTests.FEASIBLE_PROBLEM => MOI.FEASIBLE_POINT,
        MINLPTests.INFEASIBLE_PROBLEM => MOI.INFEASIBLE_POINT,
    )
    # This function tests (potentially) non-convex nonlinear programs. The tests
    # are meant to be "easy" in the sense that most NLP solvers can find the
    # same global minimum, but a test failure can sometimes be allowed.
    MINLPTests.test_nlp_expr(
        Optimizer_Uno_filtersqp;
        exclude = [
            "003_014",  # Local solution
            "004_010",  # Local solution
            "004_011",  # Local solution
            # no gradient at the initial point
            "005_010",
            "007_010",
            "008_010",  # Local solution
            # Okay to exclude forever: AmplNLWriter does not support
            # user-defined functions.
            "006_010",
        ],
        primal_target,
    )
    # This function tests convex nonlinear programs. Test failures here should
    # never be allowed, because even local NLP solvers should find the global
    # optimum.
    MINLPTests.test_nlp_cvx_expr(Optimizer_Uno_filtersqp; primal_target)
end

# This testset runs the full gamut of MOI.Test.runtests. There are a number of
# tests in here with weird edge cases, so a variety of exclusions are expected.
@testset "MathOptInterface.test" begin
    optimizer = MOI.instantiate(
        Optimizer_Uno_filtersqp;
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
            r"^test_quadratic_SecondOrderCone_basic$",
            r"^test_quadratic_nonconvex_constraint_basic$",
            r"^test_linear_DUAL_INFEASIBLE$",
            r"^test_linear_DUAL_INFEASIBLE_2$",
            r"^test_solve_TerminationStatus_DUAL_INFEASIBLE$",
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
