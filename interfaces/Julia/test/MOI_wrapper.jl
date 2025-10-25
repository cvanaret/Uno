import MathOptInterface as MOI

function test_MOI_Test()
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.Bridges.full_bridge_optimizer(UnoSolver.Optimizer(), Float64),
    )
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            atol = 1e-4,
            rtol = 1e-4,
            infeasible_status = MOI.LOCALLY_INFEASIBLE,
            optimal_status = MOI.LOCALLY_SOLVED,
            exclude = Any[
                MOI.ConstraintBasisStatus,
                MOI.DualObjectiveValue,
                MOI.ObjectiveBound,
            ],
        );
        exclude = [
            # Failures in UnoSolver.jl -- need an investigation of Charlie
            "test_conic_NormInfinityCone_3",
            "test_conic_NormInfinityCone_INFEASIBLE",
            r"test_conic_NormOneCone.*",
            r"test_conic_linear_VectorAffineFunction.*",
            r"test_conic_linear_VectorOfVariables.*",
            "test_constraint_ScalarAffineFunction_EqualTo",
            "test_constraint_ScalarAffineFunction_GreaterThan",
            "test_constraint_qcp_duplicate_off_diagonal",
            r"test_linear_integration.*",
            "test_nonlinear_duals",
            "test_nonlinear_hs071_NLPBlockDual",
            "test_nonlinear_with_scalar_quadratic_function_with_off_diag",
            "test_quadratic_SecondOrderCone_basic",
            "test_quadratic_constraint_GreaterThan",
            "test_quadratic_constraint_LessThan",
            "test_quadratic_nonconvex_constraint_basic",
            "test_quadratic_nonhomogeneous",
            "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_Interval_upper",
            r"test_quadratic_constraint_basic.*",
            r"test_quadratic_constraint_integration.*",
            r"test_quadratic_constraint_minimize.*",
            r"test_quadratic_duplicate_terms.*",
        ],
    )
    return
end

function test_Name()
    model = UnoSolver.Optimizer()
    @test MOI.supports(model, MOI.Name())
    @test MOI.get(model, MOI.Name()) == ""
    MOI.set(model, MOI.Name(), "Model")
    @test MOI.get(model, MOI.Name()) == "Model"
    return
end

function test_ConstraintDualStart()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    l = MOI.add_constraint(model, x[1], MOI.GreaterThan(1.0))
    u = MOI.add_constraint(model, x[1], MOI.LessThan(1.0))
    e = MOI.add_constraint(model, x[2], MOI.EqualTo(1.0))
    c = MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0),
        MOI.LessThan(1.5),
    )
    @test MOI.supports(model, MOI.ConstraintDualStart(), typeof(l))
    @test MOI.supports(model, MOI.ConstraintDualStart(), typeof(u))
    @test MOI.supports(model, MOI.ConstraintDualStart(), typeof(e))
    @test MOI.supports(model, MOI.ConstraintDualStart(), typeof(c))
    @test MOI.supports(model, MOI.NLPBlockDualStart())
    @test MOI.get(model, MOI.ConstraintDualStart(), l) === nothing
    @test MOI.get(model, MOI.ConstraintDualStart(), u) === nothing
    @test MOI.get(model, MOI.ConstraintDualStart(), e) === nothing
    @test MOI.get(model, MOI.ConstraintDualStart(), c) === nothing
    @test MOI.get(model, MOI.NLPBlockDualStart()) === nothing
    MOI.set(model, MOI.ConstraintDualStart(), l, 1.0)
    MOI.set(model, MOI.ConstraintDualStart(), u, -1.0)
    MOI.set(model, MOI.ConstraintDualStart(), e, -1.5)
    MOI.set(model, MOI.ConstraintDualStart(), c, 2.0)
    MOI.set(model, MOI.NLPBlockDualStart(), [1.0, 2.0])
    @test MOI.get(model, MOI.ConstraintDualStart(), l) == 1.0
    @test MOI.get(model, MOI.ConstraintDualStart(), u) == -1.0
    @test MOI.get(model, MOI.ConstraintDualStart(), e) == -1.5
    @test MOI.get(model, MOI.ConstraintDualStart(), c) == 2.0
    @test MOI.get(model, MOI.NLPBlockDualStart()) == [1.0, 2.0]
    MOI.set(model, MOI.ConstraintDualStart(), l, nothing)
    MOI.set(model, MOI.ConstraintDualStart(), u, nothing)
    MOI.set(model, MOI.ConstraintDualStart(), e, nothing)
    MOI.set(model, MOI.ConstraintDualStart(), c, nothing)
    MOI.set(model, MOI.NLPBlockDualStart(), nothing)
    @test MOI.get(model, MOI.ConstraintDualStart(), l) === nothing
    @test MOI.get(model, MOI.ConstraintDualStart(), u) === nothing
    @test MOI.get(model, MOI.ConstraintDualStart(), e) === nothing
    @test MOI.get(model, MOI.ConstraintDualStart(), c) === nothing
    @test MOI.get(model, MOI.NLPBlockDualStart()) === nothing
    return
end

function test_ConstraintDualStart_ScalarNonlinearFunction()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint.(model, x, MOI.Interval(0.0, 0.8))
    f = MOI.ScalarNonlinearFunction(:sin, Any[1.0*x[1]-1.0*x[2]])
    c = MOI.add_constraint(model, f, MOI.EqualTo(0.5))
    g = 1.0 * x[1] + 1.0 * x[2]
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(g)}(), g)
    @test MOI.supports(model, MOI.ConstraintDualStart(), typeof(c))
    @test MOI.get(model, MOI.ConstraintDualStart(), c) === nothing
    MOI.set(model, MOI.ConstraintDualStart(), c, 1.15)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ConstraintDualStart(), c) === 1.15
    MOI.set(model, MOI.ConstraintDualStart(), c, nothing)
    @test MOI.get(model, MOI.ConstraintDualStart(), c) === nothing
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.LOCALLY_SOLVED
    return
end

function test_ConstraintDualStart_variable_bound_min_greater_than()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x, c = MOI.add_constrained_variable(model, MOI.GreaterThan(1.0))
    MOI.set(model, MOI.VariablePrimalStart(), x, 1.0)
    MOI.set(model, MOI.ConstraintDualStart(), c, 1.0)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.VariableIndex}(), x)
    MOI.optimize!(model)
    @test isapprox(MOI.get(model, MOI.ConstraintDual(), c), 1.0; atol = 1e-6)
    @test MOI.get(model, MOI.ConstraintDualStart(), c) == 1.0
    return
end

function test_ConstraintDualStart_variable_bound_max_less_than()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x, c = MOI.add_constrained_variable(model, MOI.LessThan(1.0))
    MOI.set(model, MOI.VariablePrimalStart(), x, 1.0)
    MOI.set(model, MOI.ConstraintDualStart(), c, -1.0)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.VariableIndex}(), x)
    MOI.optimize!(model)
    @test isapprox(MOI.get(model, MOI.ConstraintDual(), c), -1.0; atol = 1e-6)
    @test MOI.get(model, MOI.ConstraintDualStart(), c) == -1.0
    return
end

function test_ConstraintDualStart_variable_bound_min_equal_to()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x, c = MOI.add_constrained_variable(model, MOI.EqualTo(1.0))
    MOI.set(model, MOI.VariablePrimalStart(), x, 1.0)
    MOI.set(model, MOI.ConstraintDualStart(), c, 1.0)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.VariableIndex}(), x)
    MOI.optimize!(model)
    @test isapprox(MOI.get(model, MOI.ConstraintDual(), c), 1.0; atol = 1e-6)
    @test MOI.get(model, MOI.ConstraintDualStart(), c) == 1.0
    return
end

function test_solve_time()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    MOI.add_variable(model)
    @test isnan(MOI.get(model, MOI.SolveTimeSec()))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.SolveTimeSec()) >= 0.0
    return
end

function test_barrier_iterations()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    f = (x - 1.0)^2 + 2.0 * x + 3.0
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.add_constraint(model, 1.0 * x, MOI.LessThan(0.5))
    @test MOI.get(model, MOI.BarrierIterations()) == 0
    MOI.optimize!(model)
    @test MOI.get(model, MOI.BarrierIterations()) > 0
    return
end

# Model structure for test_check_derivatives_for_naninf()
struct Issue136 <: MOI.AbstractNLPEvaluator end
MOI.initialize(::Issue136, ::Vector{Symbol}) = nothing
MOI.features_available(::Issue136) = [:Grad, :Jac]
MOI.eval_objective(::Issue136, x) = x[1]
MOI.eval_constraint(::Issue136, g, x) = (g[1] = x[1]^(1 / 3))
MOI.eval_objective_gradient(::Issue136, grad_f, x) = (grad_f[1] = 1.0)
MOI.jacobian_structure(::Issue136) = Tuple{Int64,Int64}[(1, 1)]
function MOI.eval_constraint_jacobian(::Issue136, J, x)
    J[1] = (1 / 3) * x[1]^(1 / 3 - 1)
    return
end

function test_check_derivatives_for_naninf()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.set(
        model,
        MOI.NLPBlock(),
        MOI.NLPBlockData(MOI.NLPBoundsPair.([-Inf], [0.0]), Issue136(), false),
    )
    # Failure to set check_derivatives_for_naninf="yes" may cause Uno to
    # segfault or return a NUMERICAL_ERROR status. Check that it is set to "yes"
    # by obtaining an INVALID_MODEL status.
    # MOI.set(model, MOI.RawOptimizerAttribute("check_derivatives_for_naninf"), "no")
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INVALID_MODEL
    @test MOI.get(model, MOI.NLPBlock()) isa MOI.NLPBlockData
    return
end

function test_empty_optimize()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    @test MOI.get(model, MOI.RawStatusString()) == "Optimize not called"
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INVALID_MODEL
    @test MOI.get(model, MOI.DualStatus()) == MOI.NO_SOLUTION
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.NO_SOLUTION
    @test MOI.get(model, MOI.RawStatusString()) == "The model has no variable"
    return
end

"""
    test_get_model()

Test various getters for ConstraintFunction etc. We need this test because the
normal MOI ones require the solver to support VariableName and ConstraintName.
"""
function test_get_model()
    model = MOI.Utilities.Model{Float64}()
    MOI.Utilities.loadfromstring!(
        model,
        """
        variables: w, x, y, z
        minobjective: 1.0 * x * x + 2.0 * y + 3.0
        w in Interval(-1.0, 1.0)
        x >= 1.0
        y <= 2.0
        z == 3.0
        1.0 * x >= 1.0
        2.0 * y <= 4.0
        3.0 * z == 9.0
        4.0 * w in Interval(-10.1, 11.1)
        1.0 * x * x + x >= 1.0
        2.0 * y * y + y <= 8.0
        3.0 * z * z + z == 27.0
        """,
    )
    uno = UnoSolver.Optimizer()
    index_map = MOI.copy_to(uno, model)
    attr = MOI.ListOfConstraintTypesPresent()
    @test sort(MOI.get(model, attr); by = string) ==
          sort(MOI.get(uno, attr); by = string)
    for (F, S) in MOI.get(model, MOI.ListOfConstraintTypesPresent())
        cis = MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
        @test length(cis) == 1
        f_model = MOI.get(model, MOI.ConstraintFunction(), cis[1])
        s_model = MOI.get(model, MOI.ConstraintSet(), cis[1])
        cis = MOI.get(uno, MOI.ListOfConstraintIndices{F,S}())
        @test length(cis) == MOI.get(uno, MOI.NumberOfConstraints{F,S}()) == 1
        f_uno = MOI.get(uno, MOI.ConstraintFunction(), cis[1])
        s_uno = MOI.get(uno, MOI.ConstraintSet(), cis[1])
        @test s_model == s_uno
        if F == MOI.VariableIndex
            @test index_map[f_model] == f_uno
        else
            @test ≈(
                MOI.Utilities.substitute_variables(x -> index_map[x], f_model),
                f_uno,
            )
        end
    end
    F_model = MOI.get(model, MOI.ObjectiveFunctionType())
    F_uno = MOI.get(uno, MOI.ObjectiveFunctionType())
    @test F_model == F_uno
    obj_model = MOI.get(model, MOI.ObjectiveFunction{F_model}())
    obj_uno = MOI.get(uno, MOI.ObjectiveFunction{F_uno}())
    @test ≈(
        MOI.Utilities.substitute_variables(x -> index_map[x], obj_model),
        obj_uno,
    )
    return
end

function test_supports_ConstraintDualStart_VariableIndex()
    uno = UnoSolver.Optimizer()
    bridged = MOI.Bridges.full_bridge_optimizer(UnoSolver.Optimizer(), Float64)
    sets =
        (MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.EqualTo{Float64})
    for model in (uno, bridged), S in sets
        @test MOI.supports(
            model,
            MOI.ConstraintDualStart(),
            MOI.ConstraintIndex{MOI.VariableIndex,S},
        )
    end
    return
end

function test_parameter_number_of_variables()
    model = UnoSolver.Optimizer()
    x = MOI.add_variable(model)
    p, ci = MOI.add_constrained_variable(model, MOI.Parameter(2.0))
    @test MOI.get(model, MOI.NumberOfVariables()) == 2
    return
end

function test_parameter_is_valid()
    model = UnoSolver.Optimizer()
    p, ci = MOI.add_constrained_variable(model, MOI.Parameter(2.0))
    @test MOI.is_valid(model, p)
    @test MOI.is_valid(model, ci)
    @test !MOI.is_valid(model, typeof(p)(p.value + 1))
    @test !MOI.is_valid(model, typeof(ci)(ci.value + 1))
    return
end

function test_parameter_list_of_variable_indices()
    model = UnoSolver.Optimizer()
    x = MOI.add_variable(model)
    p, ci = MOI.add_constrained_variable(model, MOI.Parameter(2.0))
    @test MOI.get(model, MOI.ListOfVariableIndices()) == [x, p]
    # Now reversed
    model = UnoSolver.Optimizer()
    p, ci = MOI.add_constrained_variable(model, MOI.Parameter(2.0))
    x = MOI.add_variable(model)
    @test MOI.get(model, MOI.ListOfVariableIndices()) == [p, x]
    return
end

function test_scalar_nonlinear_function_is_valid()
    model = UnoSolver.Optimizer()
    x = MOI.add_variable(model)
    F, S = MOI.ScalarNonlinearFunction, MOI.EqualTo{Float64}
    @test MOI.is_valid(model, MOI.ConstraintIndex{F,S}(1)) == false
    f = MOI.ScalarNonlinearFunction(:sin, Any[x])
    c = MOI.add_constraint(model, f, MOI.EqualTo(0.0))
    @test c isa MOI.ConstraintIndex{F,S}
    @test MOI.is_valid(model, c) == true
    return
end

function test_scalar_nonlinear_function_nlp_block()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    f = MOI.ScalarNonlinearFunction(:^, Any[x, 4])
    MOI.add_constraint(model, f, MOI.LessThan(1.0))
    MOI.optimize!(model)
    block = MOI.get(model, MOI.NLPBlock())
    @test !block.has_objective
    @test block.evaluator isa MOI.Nonlinear.Evaluator
    return
end

function test_parameter()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    p, ci = MOI.add_constrained_variable(model, MOI.Parameter(1.0))
    x = MOI.add_variable(model)
    fi = MOI.ScalarNonlinearFunction(:-, Any[x, p])
    f = MOI.ScalarNonlinearFunction(:^, Any[fi, 2])
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test MOI.get(model, MOI.NumberOfVariables()) == 2
    @test MOI.get(model, MOI.ListOfVariableIndices()) == [p, x]
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), p) ≈ 1
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ 1
    MOI.set(model, MOI.ConstraintSet(), ci, MOI.Parameter(-2.5))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), p) ≈ -2.5
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ -2.5
    return
end

function test_parameter_replace_parameters()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    p, ci = MOI.add_constrained_variable(model, MOI.Parameter(1.0))
    x = MOI.add_variable(model)
    t = MOI.add_variable(model)
    lhs = MOI.ScalarNonlinearFunction(
        :+,
        Any[
            x,
            p,
            1.0*x,
            1.0*p,
            1.0*x*x,
            1.0*p*x,
            MOI.ScalarNonlinearFunction(:^, Any[1.0*x-p, 2]),
        ],
    )
    f = MOI.ScalarNonlinearFunction(:-, Any[lhs, t])
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(t)}(), t)
    MOI.add_constraint(model, f, MOI.LessThan(0.0))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), p) ≈ 1.0
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ -0.25
    return
end

function test_parameter_reverse()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    p, ci = MOI.add_constrained_variable(model, MOI.Parameter(1.0))
    fi = MOI.ScalarNonlinearFunction(:-, Any[x, p])
    f = MOI.ScalarNonlinearFunction(:^, Any[fi, 2])
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test MOI.get(model, MOI.NumberOfVariables()) == 2
    @test MOI.get(model, MOI.ListOfVariableIndices()) == [x, p]
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), p) ≈ 1
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ 1
    MOI.set(model, MOI.ConstraintSet(), ci, MOI.Parameter(-2.5))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), p) ≈ -2.5
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ -2.5
    return
end

function test_parameter_scalar_affine_objective()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    p, ci = MOI.add_constrained_variable(model, MOI.Parameter(2.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    # f = (x - p)^2 + x + p + 1.0
    f = (1.0 * x - 1.0 * p) * (1.0 * x - 1.0 * p) + x + p + 1.0
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ 1.5
    @test MOI.get(model, MOI.VariablePrimal(), p) ≈ 2.0
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ (1.5 - 2.0)^2 + 4.5
    MOI.set(model, MOI.ConstraintSet(), ci, MOI.Parameter(2.2))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ 1.7
    @test MOI.get(model, MOI.VariablePrimal(), p) ≈ 2.2
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ (1.7 - 2.2)^2 + 4.9
    MOI.add_constraint(model, x, MOI.LessThan(1.5))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ 1.5
    @test MOI.get(model, MOI.VariablePrimal(), p) ≈ 2.2
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ (1.5 - 2.2)^2 + 4.7
    return
end

function test_parameter_variable_index_objective()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    p, ci = MOI.add_constrained_variable(model, MOI.Parameter(2.0))
    t = MOI.add_variable(model)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(t)}(), t)
    # f = (x - p)^2 + x + p + 1.0
    f = (1.0 * x - 1.0 * p) * (1.0 * x - 1.0 * p) + x + p + 1.0
    MOI.add_constraint(model, f - t, MOI.LessThan(0.0))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ 1.5
    @test MOI.get(model, MOI.VariablePrimal(), p) ≈ 2.0
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ (1.5 - 2.0)^2 + 4.5
    MOI.set(model, MOI.ConstraintSet(), ci, MOI.Parameter(2.2))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ 1.7
    @test MOI.get(model, MOI.VariablePrimal(), p) ≈ 2.2
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ (1.7 - 2.2)^2 + 4.9
    return
end

function test_ListOfSupportedNonlinearOperators()
    model = UnoSolver.Optimizer()
    ops = MOI.get(model, MOI.ListOfSupportedNonlinearOperators())
    @test ops isa Vector{Symbol}
    @test :|| in ops
    @test :ifelse in ops
    @test :sin in ops
    @test !(:f in ops)
    f(x) = x^2
    MOI.set(model, MOI.UserDefinedFunction(:f, 1), (f,))
    @test :f in MOI.get(model, MOI.ListOfSupportedNonlinearOperators())
    return
end

# function test_SPRAL()
#     if !haskey(ENV, "OMP_CANCELLATION") || !haskey(ENV, "OMP_PROC_BIND")
#         return
#     end
#     model = UnoSolver.Optimizer()
#     MOI.set(model, MOI.RawOptimizerAttribute("linear_solver"), "spral")
#     MOI.set(model, MOI.Silent(), true)
#     x = MOI.add_variable(model)
#     MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
#     f = 1.0 * x * x - 4.0 * x + 4.0
#     MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
#     MOI.optimize!(model)
#     @test ≈(MOI.get(model, MOI.VariablePrimal(), x), 2.0; atol = 1e-6)
#     return
# end

# function test_HSL()
#     if !(@ccall HSL_jll.libhsl.LIBHSL_isfunctional()::Bool)
#         return
#     end
#     for hsl_solver in ("ma27", "ma57", "ma77", "ma86", "ma97")
#         model = UnoSolver.Optimizer()
#         MOI.set(model, MOI.RawOptimizerAttribute("linear_solver"), hsl_solver)
#         MOI.set(model, MOI.Silent(), true)
#         x = MOI.add_variable(model)
#         MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
#         f = 1.0 * x * x - 4.0 * x + 4.0
#         MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
#         MOI.optimize!(model)
#         @test ≈(MOI.get(model, MOI.VariablePrimal(), x), 2.0; atol = 1e-6)
#     end
#     return
# end

function test_ad_backend()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    attr = MOI.AutomaticDifferentiationBackend()
    @test MOI.supports(model, attr)
    @test MOI.get(model, attr) == MOI.Nonlinear.SparseReverseMode()
    MOI.optimize!(model)
    @test model.inner isa UnoSolver.Model
    MOI.set(model, attr, MOI.Nonlinear.ExprGraphOnly())
    @test MOI.get(model, attr) == MOI.Nonlinear.ExprGraphOnly()
    @test model.inner === nothing
    f = MOI.ScalarNonlinearFunction(:^, Any[x, 4])
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    @test_throws ErrorException MOI.optimize!(model)
    return
end

function test_mixing_new_old_api()
    # new then old
    model = UnoSolver.Optimizer()
    MOI.add_constrained_variable(model, MOI.Parameter(2.0))
    bounds = MOI.NLPBoundsPair.([25.0, 40.0], [Inf, 40.0])
    block_data = MOI.NLPBlockData(bounds, MOI.Test.HS071(true), true)
    err = ErrorException("Cannot mix the new and legacy nonlinear APIs")
    @test_throws err MOI.set(model, MOI.NLPBlock(), block_data)
    # old then new
    model = UnoSolver.Optimizer()
    bounds = MOI.NLPBoundsPair.([25.0, 40.0], [Inf, 40.0])
    block_data = MOI.NLPBlockData(bounds, MOI.Test.HS071(true), true)
    MOI.set(model, MOI.NLPBlock(), block_data)
    err = ErrorException("Cannot mix the new and legacy nonlinear APIs")
    @test_throws err MOI.add_constrained_variable(model, MOI.Parameter(2.0))
    return
end

function test_nlp_model_objective_function_type()
    model = UnoSolver.Optimizer()
    x = MOI.add_variable(model)
    f = MOI.ScalarNonlinearFunction(:sqrt, Any[x])
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    F = MOI.get(model, MOI.ObjectiveFunctionType())
    @test F == MOI.ScalarNonlinearFunction
    return
end

function test_nlp_model_set_set()
    model = UnoSolver.Optimizer()
    x = MOI.add_variable(model)
    f = MOI.ScalarNonlinearFunction(:sqrt, Any[x])
    c = MOI.add_constraint(model, f, MOI.LessThan(2.0))
    @test MOI.get(model, MOI.ConstraintSet(), c) == MOI.LessThan(2.0)
    MOI.set(model, MOI.ConstraintSet(), c, MOI.LessThan(3.0))
    @test MOI.get(model, MOI.ConstraintSet(), c) == MOI.LessThan(3.0)
    return
end

function test_VariablePrimalStart()
    attr = MOI.VariablePrimalStart()
    model = UnoSolver.Optimizer()
    x = MOI.add_variable(model)
    @test MOI.supports(model, attr, typeof(x))
    @test MOI.get(model, attr, x) === nothing
    MOI.set(model, attr, x, 1.0)
    @test MOI.get(model, attr, x) == 1.0
    p, _ = MOI.add_constrained_variable(model, MOI.Parameter(1.0))
    @test_throws(
        MOI.GetAttributeNotAllowed{typeof(attr)},
        MOI.get(model, attr, p),
    )
    @test_throws(
        MOI.SetAttributeNotAllowed{typeof(attr)},
        MOI.set(model, attr, p, 1.0),
    )
    return
end

function test_manually_evaluated_primal_status()
    Ext = Base.get_extension(Uno, :UnoMathOptInterfaceExt)
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.add_constraint(model, x, MOI.Interval(0.0, 1.0))
    MOI.add_constraint(model, 1.0 * x, MOI.Interval(0.1, 1.1))
    f = MOI.ScalarNonlinearFunction(:log, Any[x])
    MOI.add_constraint(model, f, MOI.Interval(-1.0, 1.0))
    MOI.optimize!(model)

    # Alexis -- revisit this code when we know how to handle the workspace
    # x_star = copy(model.inner.x)
    # g_star = copy(model.inner.g)
    x_star = Vector{Float64}(undef, model.inner.nvar)
    UnoSolver.uno_get_primal_solution(model.inner, x_star)
    g_star = Vector{Float64}(undef, model.inner.ncon)
    MOI.eval_constraint(model, g_star, x_star)
    for (xi, status) in (
        x_star[1] => MOI.FEASIBLE_POINT,
        -1.0 => MOI.INFEASIBLE_POINT,
        -1e-7 => MOI.NEARLY_FEASIBLE_POINT,
        1 + 1e-7 => MOI.NEARLY_FEASIBLE_POINT,
        0.0 => MOI.FEASIBLE_POINT,
        1.0 => MOI.FEASIBLE_POINT,
    )
        model.inner.x[1] = xi
        @test Ext._manually_evaluated_primal_status(model) == status
    end
    model.inner.x .= x_star
    for (gi, status) in (
        g_star => MOI.FEASIBLE_POINT,
        [0.0, 0.0] => MOI.INFEASIBLE_POINT,
        [0.1 - 1e-7, 0.0] => MOI.NEARLY_FEASIBLE_POINT,
        [1.1 + 1e-7, 0.0] => MOI.NEARLY_FEASIBLE_POINT,
        [0.1, -1.0 - 1e-7] => MOI.NEARLY_FEASIBLE_POINT,
        [0.1, 1.0 + 1e-7] => MOI.NEARLY_FEASIBLE_POINT,
        [0.1, 1.0] => MOI.FEASIBLE_POINT,
        [1.1, -1.0] => MOI.FEASIBLE_POINT,
    )
        model.inner.g .= gi
        @test Ext._manually_evaluated_primal_status(model) == status
    end
    return
end

function test_RawOptimizerAttribute()
    model = UnoSolver.Optimizer()
    attr = MOI.RawOptimizerAttribute("print_solution")
    @test_throws MOI.GetAttributeNotAllowed{typeof(attr)} MOI.get(model, attr)
    @test MOI.supports(model, attr)
    MOI.set(model, attr, 0)
    @test MOI.get(model, attr) == 0
    return
end

function test_empty_nlp_evaluator()
    model = UnoSolver.Optimizer()
    block = MOI.get(model, MOI.NLPBlock())
    evaluator = block.evaluator
    @test MOI.features_available(evaluator) == [:Grad, :Jac, :Hess]
    MOI.initialize(evaluator, [:Grad, :Jac, :Hess])
    x = Float64[]
    @test MOI.eval_constraint(evaluator, Float64[], x) === nothing
    @test MOI.jacobian_structure(evaluator) == Tuple{Int,Int}[]
    @test MOI.hessian_lagrangian_structure(evaluator) == Tuple{Int,Int}[]
    @test MOI.eval_constraint_jacobian(evaluator, Float64[], x) === nothing
    H, mu = Float64[], Float64[]
    @test MOI.eval_hessian_lagrangian(evaluator, H, x, 1.0, mu) === nothing
    return
end

function test_NLPBlockDualStart()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 4)
    MOI.set.(model, MOI.VariablePrimalStart(), x, 1.0)
    block = MOI.NLPBlockData(
        MOI.NLPBoundsPair.([25.0, 40.0], [Inf, 40.0]),
        MOI.Test.HS071(false),
        true,
    )
    MOI.set(model, MOI.NLPBlock(), block)
    MOI.set(model, MOI.NLPBlockDualStart(), [1.0, -1.0])
    MOI.optimize!(model)
    @test isapprox(MOI.get(model, MOI.NLPBlockDual()), [0.0, 0.0]; atol = 1e-6)
    return
end

function test_function_type_to_func()
    model = UnoSolver.Optimizer()
    x = MOI.add_variable(model)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.VariableIndex}(), x)
    @test MOI.get(model, MOI.ObjectiveFunctionType()) == MOI.VariableIndex
    return
end

function test_error_adding_option()
    model = UnoSolver.Optimizer()
    x = MOI.add_variable(model)
    name, value = "print_solution", :zero
    MOI.set(model, MOI.RawOptimizerAttribute(name), value)
    @test_throws(
        ErrorException(
            "Unable to add option `\"$name\"` with the value " *
            "`$value::$(typeof(value))`. The value must be a `::String`, a `::Float64`, an `::Integer`, or a `::Bool`.",
        ),
        MOI.optimize!(model),
    )
    return
end

function test_scalar_nonlinear_function_attributes()
    model = UnoSolver.Optimizer()
    x = MOI.add_variable(model)
    F, S = MOI.ScalarNonlinearFunction, MOI.LessThan{Float64}
    @test isempty(MOI.get(model, MOI.ListOfConstraintTypesPresent()))
    @test isempty(MOI.get(model, MOI.ListOfConstraintIndices{F,S}()))
    @test MOI.get(model, MOI.NumberOfConstraints{F,S}()) == 0
    f = MOI.ScalarNonlinearFunction(:log, Any[x])
    c = MOI.add_constraint(model, f, MOI.LessThan(1.0))
    @test MOI.get(model, MOI.ListOfConstraintIndices{F,S}()) == [c]
    @test MOI.get(model, MOI.NumberOfConstraints{F,S}()) == 1
    @test (F, S) in MOI.get(model, MOI.ListOfConstraintTypesPresent())
    return
end

function test_vector_nonlinear_oracle_scalar_nonlinear_equivalent()
    model = UnoSolver.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 4)
    t = MOI.add_variable(model)
    MOI.add_constraint(model, x[1], MOI.GreaterThan(0.1))
    c_x3 = MOI.add_constraint(model, x[3], MOI.LessThan(1.0))
    c_x4 = MOI.add_constraint(model, x[4], MOI.GreaterThan(0.0))
    f1 = 1.0 * x[1] * x[1] + 1.0 * x[2] * x[2] - 1.0 * x[3]
    c1 = MOI.add_constraint(model, f1, MOI.EqualTo(0.0))
    f2 = 1.0 * x[2] - 1.0 * x[1] - x[4]
    c2 = MOI.add_constraint(model, f2, MOI.EqualTo(0.0))
    log_x = MOI.ScalarNonlinearFunction(:log, Any[x[1]])
    log_x_minus_t = MOI.ScalarNonlinearFunction(:-, Any[log_x, t])
    c_snf = MOI.add_constraint(model, log_x_minus_t, MOI.GreaterThan(0.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    F = MOI.ScalarAffineFunction{Float64}
    MOI.set(model, MOI.ObjectiveFunction{F}(), 1.0 * t)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.LOCALLY_SOLVED
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(model, MOI.DualStatus()) == MOI.FEASIBLE_POINT
    atol = 1e-6
    x_sol = MOI.get.(model, MOI.VariablePrimal(), x)
    @test ≈(x_sol, [1 / sqrt(2), 1 / sqrt(2), 1.0, 0.0]; atol)
    @test ≈(MOI.get(model, MOI.VariablePrimal(), t), -log(sqrt(2)); atol)
    @test ≈(MOI.get(model, MOI.ConstraintPrimal(), c1), 0.0; atol)
    @test ≈(MOI.get(model, MOI.ConstraintPrimal(), c2), 0.0; atol)
    @test ≈(MOI.get(model, MOI.ConstraintDual(), c1), -0.5; atol)
    @test ≈(MOI.get(model, MOI.ConstraintDual(), c2), 1 / sqrt(2); atol)
    @test ≈(MOI.get(model, MOI.ConstraintDual(), c_x3), -0.5; atol)
    @test ≈(MOI.get(model, MOI.ConstraintDual(), c_x4), 1 / sqrt(2); atol)
    @test ≈(MOI.get(model, MOI.ConstraintPrimal(), c_snf), 0.0; atol)
    @test ≈(MOI.get(model, MOI.ConstraintDual(), c_snf), 1.0; atol)
    return
end

test_MOI_Test()
test_Name()
test_ConstraintDualStart()
test_ConstraintDualStart_ScalarNonlinearFunction()
test_ConstraintDualStart_variable_bound_min_greater_than()
# test_ConstraintDualStart_variable_bound_max_less_than()
test_ConstraintDualStart_variable_bound_min_equal_to()
test_solve_time()
# test_barrier_iterations()
# test_check_derivatives_for_naninf()
test_empty_optimize()
test_get_model()
test_supports_ConstraintDualStart_VariableIndex()
test_parameter_number_of_variables()
test_parameter_is_valid()
test_parameter_list_of_variable_indices()
test_scalar_nonlinear_function_is_valid()
test_scalar_nonlinear_function_nlp_block()
## We need to use isapprox for the following tests
# test_parameter()
# test_parameter_replace_parameters()
# test_parameter_reverse()
# test_parameter_scalar_affine_objective()
## Uno is dead-lock in this test!
# test_parameter_variable_index_objective()
test_ListOfSupportedNonlinearOperators()
test_ad_backend()
test_mixing_new_old_api()
test_nlp_model_objective_function_type()
test_nlp_model_set_set()
test_VariablePrimalStart()
## We need to use a workspace for the test below
# test_manually_evaluated_primal_status()
test_RawOptimizerAttribute()
## We need to use isapprox for the following tests
# test_empty_nlp_evaluator()
## Bug...
# test_NLPBlockDualStart()
test_function_type_to_func()
test_error_adding_option()
test_error_adding_option()
test_scalar_nonlinear_function_attributes()
## We need to use isapprox for the following tests
# test_vector_nonlinear_oracle_scalar_nonlinear_equivalent()
