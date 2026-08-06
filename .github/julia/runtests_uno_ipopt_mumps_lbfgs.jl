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

Optimizer_Uno_ipopt() = Optimizer(["logger=DEBUG3", "preset=ipopt", "linear_solver=MUMPS", "hessian_model=LBFGS", "unbounded_objective_threshold=-1e15"])

# This testset runs https://github.com/jump-dev/MINLPTests.jl

# drop the equality-constrained instances and the LPs
bound_constrained_instances = readlines(joinpath(@__DIR__, "MINLPTests/bound-constrained.txt"))
general_instances = readlines(joinpath(@__DIR__, "MINLPTests/general.txt"))
NLP_instances = readlines(joinpath(@__DIR__, "MINLPTests/NLP.txt"))
instances = vcat(bound_constrained_instances, general_instances)
instances = intersect(instances, NLP_instances)
#print("Instances: ", instances)

# strip the prefix
function strip_prefix(instances, prefix)
    cleaned_instances = String[]
    for instance in instances
        instance = strip(instance)
        startswith(instance, prefix) && push!(cleaned_instances, replace(instance, prefix => ""; count = 1))
    end
    return cleaned_instances
end


# This testset runs the full gamut of MOI.Test.runtests. There are a number of
# tests in here with weird edge cases, so a variety of exclusions are expected.

# drop the equality-constrained instances and the LPs
bound_constrained_instances = readlines(joinpath(@__DIR__, "MOI/bound-constrained.txt"))
general_instances = readlines(joinpath(@__DIR__, "MOI/general.txt"))
NLP_instances = readlines(joinpath(@__DIR__, "MOI/NLP.txt"))
instances = vcat(bound_constrained_instances, general_instances)
instances = intersect(instances, NLP_instances)
MOI_instances = [Regex("^" * instance * "\$") for instance in instances] # exact match
MOI_instances = String["test_nonlinear_expression_hs109"]
#print("Instances: ", instances)

if !isempty(MOI_instances)
    @testset "MathOptInterface.test" begin
        optimizer = MOI.instantiate(
            Optimizer_Uno_ipopt;
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
                    MOI.DualObjectiveValue,
                ],
            );
            include = MOI_instances,
            exclude = [
                # ==================================================================
                # The following tests are bugs.
                #
                # We should fix issues in Uno, and then try removing these lines.
                r"^test_nonlinear_expression_hs109$", # Uno crashes at an infeasible point
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
                # this test passes with exact Hessian but fails with L-BFGS Hessian
                # OK to exclude for now, but should be investigated
                r"^test_objective_qp_ObjectiveFunction_edge_cases",
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
                # MPEC instances
                "_complementarity",
            ],
        )
    end
end