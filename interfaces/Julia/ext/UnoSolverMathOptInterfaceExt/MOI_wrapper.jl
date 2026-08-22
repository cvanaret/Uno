# Copyright (c) 2013: Iain Dunning, Miles Lubin, and contributors
# 2025: Adapted for UnoSolver.jl by Alexis Montoison and Charlie Vanaret
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.


mutable struct _VectorNonlinearOracleCache
    set::MOI.VectorNonlinearOracle{Float64}
    x::Vector{Float64}
    start::Union{Nothing,Vector{Float64}}
    eval_f_timer::Float64
    eval_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64

    function _VectorNonlinearOracleCache(set::MOI.VectorNonlinearOracle{Float64})
        return new(set, zeros(set.input_dimension), nothing, 0.0, 0.0, 0.0)
    end
end

"""
    Optimizer()

Create a new Uno optimizer.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Union{Nothing,UnoSolver.Model}
    solver::Union{Nothing,UnoSolver.Solver}
    name::String
    invalid_model::Bool
    silent::Bool
    options::Dict{String,Any}
    sense::MOI.OptimizationSense
    model::MOI.Nonlinear.ModelWithQuad{Float64,MOI.Nonlinear.Model}
    variable_primal_start::Vector{Union{Nothing,Float64}}
    nlp_data::MOI.NLPBlockData
    # Whether `nlp_data` was set through the legacy `MOI.NLPBlock` API, in
    # which case it must not be rebuilt from the inner nonlinear model.
    uses_nlp_block::Bool
    nlp_dual_start::Union{Nothing,Vector{Float64}}
    mult_g_nlp::Dict{MOI.Nonlinear.ConstraintIndex,Float64}
    # The evaluator of `model`, rebuilt in `_setup_model`.
    evaluator::Union{Nothing,MOI.Nonlinear.EvaluatorWithQuad{Float64}}
    ad_backend::MOI.Nonlinear.AbstractAutomaticDifferentiation
    vector_nonlinear_oracle_constraints::Vector{Tuple{MOI.VectorOfVariables,_VectorNonlinearOracleCache}}
    jrows::Vector{Cint}
    jcols::Vector{Cint}
    hrows::Vector{Cint}
    hcols::Vector{Cint}
    needs_new_inner::Bool
    hess_available::Bool
    jprod_available::Bool
    jtprod_available::Bool
    hprod_available::Bool
    problem_type::String

    function Optimizer(; kwargs...)
        option_dict = Dict{String,Any}()
        for (name, value) in kwargs
            option_dict[name |> string] = value
        end
        return new(
            nothing,
            nothing,
            "",
            false,
            false,
            option_dict,
            MOI.FEASIBILITY_SENSE,
            MOI.Nonlinear.ModelWithQuad(MOI.Nonlinear.Model()),
            Union{Nothing,Float64}[],
            MOI.NLPBlockData([], _EmptyNLPEvaluator(), false),
            false,
            nothing,
            Dict{MOI.Nonlinear.ConstraintIndex,Float64}(),
            nothing,
            MOI.Nonlinear.SparseReverseMode(),
            Tuple{MOI.VectorOfVariables,_VectorNonlinearOracleCache}[],
            Cint[],
            Cint[],
            Cint[],
            Cint[],
            true,
            false,
            false,
            false,
            false,
            "",
        )
    end
end

const _SETS = Union{
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.EqualTo{Float64},
    MOI.Interval{Float64},
}

const _FUNCTIONS = Union{
    MOI.ScalarAffineFunction{Float64},
    MOI.ScalarQuadraticFunction{Float64},
    MOI.ScalarNonlinearFunction,
}

MOI.get(::Optimizer, ::MOI.SolverVersion) = UnoSolver.uno_version() |> string

### _EmptyNLPEvaluator

struct _EmptyNLPEvaluator <: MOI.AbstractNLPEvaluator end

MOI.features_available(::_EmptyNLPEvaluator) = [:Grad, :Jac, :Hess, :JacVec, :HessVec]
MOI.initialize(::_EmptyNLPEvaluator, ::Any) = nothing
MOI.eval_constraint(::_EmptyNLPEvaluator, g, x) = nothing
MOI.jacobian_structure(::_EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::_EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.eval_constraint_jacobian(::_EmptyNLPEvaluator, J, x) = nothing
MOI.eval_hessian_lagrangian(::_EmptyNLPEvaluator, H, x, σ, μ) = nothing
MOI.eval_constraint_jacobian_product(d::_EmptyNLPEvaluator, y, x, w) = nothing
MOI.eval_constraint_jacobian_transpose_product(::_EmptyNLPEvaluator, y, x, w) = nothing
MOI.eval_hessian_lagrangian_product(::_EmptyNLPEvaluator, H, x, v, σ, μ) = nothing

function MOI.empty!(model::Optimizer)
    model.inner = nothing
    model.solver = nothing
    # SKIP: model.name
    model.invalid_model = false
    # SKIP: model.silent
    # SKIP: model.options
    model.sense = MOI.FEASIBILITY_SENSE
    model.model = MOI.Nonlinear.ModelWithQuad(MOI.Nonlinear.Model())
    empty!(model.variable_primal_start)
    model.nlp_data = MOI.NLPBlockData([], _EmptyNLPEvaluator(), false)
    model.uses_nlp_block = false
    model.nlp_dual_start = nothing
    empty!(model.mult_g_nlp)
    model.evaluator = nothing
    # SKIP: model.ad_backend
    empty!(model.vector_nonlinear_oracle_constraints)
    model.jrows = Cint[]
    model.jcols = Cint[]
    model.hrows = Cint[]
    model.hcols = Cint[]
    model.needs_new_inner = true
    model.hess_available = false
    model.jprod_available = false
    model.jtprod_available = false
    model.hprod_available = false
    model.problem_type = ""
    return
end

function MOI.is_empty(model::Optimizer)
    return MOI.is_empty(model.model.variables) &&
           isempty(model.variable_primal_start) &&
           model.nlp_data.evaluator isa _EmptyNLPEvaluator &&
           model.sense == MOI.FEASIBILITY_SENSE &&
           isempty(model.vector_nonlinear_oracle_constraints)
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(model, src)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Uno"

# The nonlinear model is `model.model.inner` and always exists; this guard
# only rejects mixing it with the legacy `MOI.NLPBlock` API.
function _check_no_nlp_block(model)
    if model.uses_nlp_block
        error("Cannot mix the new and legacy nonlinear APIs")
    end
    return
end

function MOI.supports_add_constrained_variable(
    ::Optimizer,
    ::Type{MOI.Parameter{Float64}},
)
    return true
end

function MOI.add_constrained_variable(
    model::Optimizer,
    set::MOI.Parameter{Float64},
)
    model.inner = nothing
    model.solver = nothing
    _check_no_nlp_block(model)
    return MOI.add_constrained_variable(model.model, set)
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Parameter{Float64}},
)
    return MOI.is_valid(model.model, ci)
end

function MOI.get(
    model::Optimizer,
    attr::Union{
        MOI.ListOfConstraintIndices{F,S},
        MOI.NumberOfConstraints{F,S},
    },
) where {F<:MOI.VariableIndex,S<:MOI.Parameter{Float64}}
    return MOI.get(model.model, attr)
end

function MOI.set(
    model::Optimizer,
    attr::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Parameter{Float64}},
    set::MOI.Parameter{Float64},
)
    MOI.set(model.model, attr, ci, set)
    return
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:Union{MOI.VariableIndex,_FUNCTIONS}},
    ::Type{<:_SETS},
)
    return true
end

### MOI.ListOfConstraintTypesPresent

function _add_scalar_nonlinear_constraints(ret, nlp_model::MOI.Nonlinear.Model)
    for v in values(nlp_model.constraints)
        F, S = MOI.ScalarNonlinearFunction, typeof(v.set)
        if !((F, S) in ret)
            push!(ret, (F, S))
        end
    end
    return
end

function MOI.get(model::Optimizer, attr::MOI.ListOfConstraintTypesPresent)
    ret = MOI.get(model.model.variables, attr)
    append!(ret, MOI.get(model.model, attr))
    _add_scalar_nonlinear_constraints(ret, model.model.inner)
    if !isempty(model.vector_nonlinear_oracle_constraints)
        push!(ret, (MOI.VectorOfVariables, MOI.VectorNonlinearOracle{Float64}))
    end
    if !isempty(model.model.qp.parameters)
        push!(ret, (MOI.VariableIndex, MOI.Parameter{Float64}))
    end
    return ret
end

### MOI.Name

MOI.supports(::Optimizer, ::MOI.Name) = true

function MOI.set(model::Optimizer, ::MOI.Name, value::String)
    model.name = value
    return
end

MOI.get(model::Optimizer, ::MOI.Name) = model.name

### MOI.Silent

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(model::Optimizer, ::MOI.Silent, value)
    model.silent = value
    return
end

MOI.get(model::Optimizer, ::MOI.Silent) = model.silent

### MOI.TimeLimitSec

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value::Real)
    MOI.set(model, MOI.RawOptimizerAttribute("time_limit"), Float64(value))
    return
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, ::Nothing)
    delete!(model.options, "time_limit")
    return
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    return get(model.options, "time_limit", nothing)
end

### MOI.RawOptimizerAttribute

MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute) = true

function MOI.set(model::Optimizer, p::MOI.RawOptimizerAttribute, value)
    model.options[p.name] = value
    # No need to reset model.inner and model.solver, because this gets handled in optimize!.
    return
end

function MOI.get(model::Optimizer, p::MOI.RawOptimizerAttribute)
    if !haskey(model.options, p.name)
        msg = "RawOptimizerAttribute with name $(p.name) is not already set."
        throw(MOI.GetAttributeNotAllowed(p, msg))
    end
    return model.options[p.name]
end

### Variables

"""
    column(x::MOI.VariableIndex)

Return the column associated with a variable.
"""
column(x::MOI.VariableIndex) = x.value

function MOI.add_variable(model::Optimizer)
    push!(model.variable_primal_start, nothing)
    model.inner = nothing
    model.solver = nothing
    return MOI.add_variable(model.model)
end

function MOI.is_valid(model::Optimizer, x::MOI.VariableIndex)
    return MOI.is_valid(model.model, x)
end

function MOI.get(model::Optimizer, attr::MOI.ListOfVariableIndices)
    return MOI.get(model.model, attr)
end

function MOI.get(model::Optimizer, attr::MOI.NumberOfVariables)
    return MOI.get(model.model, attr)
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    return MOI.is_valid(model.model.variables, ci)
end

function MOI.get(
    model::Optimizer,
    attr::Union{
        MOI.NumberOfConstraints{MOI.VariableIndex,<:_SETS},
        MOI.ListOfConstraintIndices{MOI.VariableIndex,<:_SETS},
    },
)
    return MOI.get(model.model.variables, attr)
end

function MOI.get(
    model::Optimizer,
    attr::Union{MOI.ConstraintFunction,MOI.ConstraintSet},
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    return MOI.get(model.model.variables, attr, c)
end

function MOI.add_constraint(model::Optimizer, x::MOI.VariableIndex, set::_SETS)
    index = MOI.add_constraint(model.model.variables, x, set)
    model.inner = nothing
    model.solver = nothing
    return index
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
    set::MOI.LessThan{Float64},
)
    MOI.set(model.model.variables, MOI.ConstraintSet(), ci, set)
    if !isnothing(model.inner) && !model.needs_new_inner
        vi = ci.value
        UnoSolver.uno_set_variable_upper_bound(model.inner, vi, model.model.variables.upper[vi])
    end
    model.solver = nothing
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
    set::MOI.GreaterThan{Float64},
)
    MOI.set(model.model.variables, MOI.ConstraintSet(), ci, set)
    if !isnothing(model.inner) && !model.needs_new_inner
        vi = ci.value
        UnoSolver.uno_set_variable_lower_bound(model.inner, vi, model.model.variables.lower[vi])
    end
    model.solver = nothing
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
    set::S,
) where {
    S <: Union{
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    },
}
    MOI.set(model.model.variables, MOI.ConstraintSet(), ci, set)
    if !isnothing(model.inner) && !model.needs_new_inner
        vi = ci.value
        UnoSolver.uno_set_variable_lower_bound(model.inner, vi, model.model.variables.lower[vi])
        UnoSolver.uno_set_variable_upper_bound(model.inner, vi, model.model.variables.upper[vi])
    end
    model.solver = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    MOI.delete(model.model.variables, ci)
    model.inner = nothing
    model.solver = nothing
    return
end

### ScalarAffineFunction and ScalarQuadraticFunction constraints

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{F,<:_SETS},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    return MOI.is_valid(model.model, ci)
end

function MOI.add_constraint(
    model::Optimizer,
    func::Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
    set::_SETS,
)
    index = MOI.add_constraint(model.model, func, set)
    model.inner = nothing
    model.solver = nothing
    return index
end

function MOI.get(
    model::Optimizer,
    attr::Union{MOI.NumberOfConstraints{F,S},MOI.ListOfConstraintIndices{F,S}},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
    S<:_SETS,
}
    return MOI.get(model.model, attr)
end

function MOI.get(
    model::Optimizer,
    attr::Union{MOI.ConstraintFunction,MOI.ConstraintSet},
    c::MOI.ConstraintIndex{F,<:_SETS},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    return MOI.get(model.model, attr, c)
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{F,MOI.LessThan{Float64}},
    set::MOI.LessThan{Float64},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    MOI.set(model.model, MOI.ConstraintSet(), ci, set)
    if !isnothing(model.inner) && !model.needs_new_inner
        UnoSolver.uno_set_constraint_upper_bound(model.inner, row(model, ci), model.model.qp.g_U[ci.value])
    end
    model.solver = nothing
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{F,MOI.GreaterThan{Float64}},
    set::MOI.GreaterThan{Float64},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    MOI.set(model.model, MOI.ConstraintSet(), ci, set)
    if !isnothing(model.inner) && !model.needs_new_inner
        UnoSolver.uno_set_constraint_lower_bound(model.inner, row(model, ci), model.model.qp.g_L[ci.value])
    end
    model.solver = nothing
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{F,S},
    set::S,
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
    S <: Union{
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    },
}
    MOI.set(model.model, MOI.ConstraintSet(), ci, set)
    if !isnothing(model.inner) && !model.needs_new_inner
        UnoSolver.uno_set_constraint_lower_bound(model.inner, row(model, ci), model.model.qp.g_L[ci.value])
        UnoSolver.uno_set_constraint_upper_bound(model.inner, row(model, ci), model.model.qp.g_U[ci.value])
    end
    model.solver = nothing
    return
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ConstraintDualStart,
    ::Type{<:MOI.ConstraintIndex{F,<:_SETS}},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    return true
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDualStart,
    c::MOI.ConstraintIndex{F,<:_SETS},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    return MOI.get(model.model, attr, c)
end

function MOI.set(
    model::Optimizer,
    attr::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{F,<:_SETS},
    value::Union{Real,Nothing},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    MOI.throw_if_not_valid(model, ci)
    MOI.set(model.model, attr, ci, value)
    # No need to reset model.inner and model.solver, because this gets handled in optimize!.
    return
end

### ScalarNonlinearFunction

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,<:_SETS},
)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    return MOI.is_valid(model.model.inner, index)
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.ScalarNonlinearFunction,
    s::_SETS,
)
    _check_no_nlp_block(model)
    index = MOI.Nonlinear.add_constraint(model.model, f, s)
    model.inner = nothing
    model.solver = nothing
    return MOI.ConstraintIndex{typeof(f),typeof(s)}(index.value)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ListOfConstraintIndices{F,S},
) where {F<:MOI.ScalarNonlinearFunction,S<:_SETS}
    ret = MOI.ConstraintIndex{F,S}[]
    for (k, v) in model.model.inner.constraints
        if v.set isa S
            push!(ret, MOI.ConstraintIndex{F,S}(k.value))
        end
    end
    return ret
end

function MOI.get(
    model::Optimizer,
    attr::MOI.NumberOfConstraints{F,S},
) where {F<:MOI.ScalarNonlinearFunction,S<:_SETS}
    return count(v.set isa S for v in values(model.model.inner.constraints))
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction},
)
    return true
end

function MOI.set(
    model::Optimizer,
    attr::MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction},
    func::MOI.ScalarNonlinearFunction,
)
    _check_no_nlp_block(model)
    MOI.Nonlinear.set_objective(model.model, func)
    model.inner = nothing
    model.solver = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,<:_SETS},
)
    MOI.throw_if_not_valid(model, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    return model.model.inner[index].set
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,MOI.LessThan{Float64}},
    set::MOI.LessThan{Float64},
)
    MOI.throw_if_not_valid(model, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    func = model.model.inner[index].expression
    model.model.inner.constraints[index] = MOI.Nonlinear.Constraint(func, set)
    if !isnothing(model.inner) && !model.needs_new_inner
        UnoSolver.uno_set_constraint_upper_bound(model.inner, row(model, ci), model.nlp_data.constraint_bounds.upper[ci.value])
    end
    model.solver = nothing
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,MOI.GreaterThan{Float64}},
    set::MOI.GreaterThan{Float64},
)
    MOI.throw_if_not_valid(model, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    func = model.model.inner[index].expression
    model.model.inner.constraints[index] = MOI.Nonlinear.Constraint(func, set)
    if !isnothing(model.inner) && !model.needs_new_inner
        UnoSolver.uno_set_constraint_lower_bound(model.inner, row(model, ci), model.nlp_data.constraint_bounds.lower[ci.value])
    end
    model.solver = nothing
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
    set::S,
) where {
    S <: Union{
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    },
}
    MOI.throw_if_not_valid(model, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    func = model.model.inner[index].expression
    model.model.inner.constraints[index] = MOI.Nonlinear.Constraint(func, set)
    if !isnothing(model.inner) && !model.needs_new_inner
        UnoSolver.uno_set_constraint_lower_bound(model.inner, row(model, ci), model.nlp_data.constraint_bounds.lower[ci.value])
        UnoSolver.uno_set_constraint_upper_bound(model.inner, row(model, ci), model.nlp_data.constraint_bounds.upper[ci.value])
    end
    model.solver = nothing
    return
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ConstraintDualStart,
    ::Type{<:MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,<:_SETS}},
)
    return true
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,<:_SETS},
)
    MOI.throw_if_not_valid(model, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    return get(model.mult_g_nlp, index, nothing)
end

function MOI.set(
    model::Optimizer,
    attr::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,<:_SETS},
    value::Union{Real,Nothing},
)
    MOI.throw_if_not_valid(model, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    if value === nothing
        delete!(model.mult_g_nlp, index)
    else
        model.mult_g_nlp[index] = convert(Float64, value)
    end
    # No need to reset model.inner and model.solver, because this gets handled in optimize!.
    return
end

### MOI.VectorOfVariables in MOI.VectorNonlinearOracle{Float64}

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorOfVariables},
    ::Type{MOI.VectorNonlinearOracle{Float64}},
)
    return true
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        MOI.VectorNonlinearOracle{Float64},
    },
)
    return 1 <= ci.value <= length(model.vector_nonlinear_oracle_constraints)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ListOfConstraintIndices{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    n = length(model.vector_nonlinear_oracle_constraints)
    return MOI.ConstraintIndex{F,S}.(1:n)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.NumberOfConstraints{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    return length(model.vector_nonlinear_oracle_constraints)
end

function MOI.add_constraint(
    model::Optimizer,
    f::F,
    s::S,
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    model.inner = nothing
    model.solver = nothing
    cache = _VectorNonlinearOracleCache(s)
    push!(model.vector_nonlinear_oracle_constraints, (f, cache))
    n = length(model.vector_nonlinear_oracle_constraints)
    return MOI.ConstraintIndex{F,S}(n)
end

function row(
    model::Optimizer,
    ci::MOI.ConstraintIndex{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    offset = length(model.model)
    for i in 1:(ci.value-1)
        _, s = model.vector_nonlinear_oracle_constraints[i]
        offset += s.set.output_dimension
    end
    _, s = model.vector_nonlinear_oracle_constraints[ci.value]
    return offset .+ (1:s.set.output_dimension)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    f, _ = model.vector_nonlinear_oracle_constraints[ci.value]
    return MOI.get.(model, MOI.VariablePrimal(attr.result_index), f.variables)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    f, s = model.vector_nonlinear_oracle_constraints[ci.value]
    J = zeros(length(s.set.jacobian_structure))
    x = [MOI.get(model, MOI.VariablePrimal(), xi) for xi in f.variables]
    s.set.eval_jacobian(J, x)
    λ = Float64[
        UnoSolver.uno_get_constraint_dual_solution_component(model.solver, r)
        for r in row(model, ci)
    ]
    dual = zeros(MOI.dimension(s.set))
    sign = _dual_multiplier(model)
    for ((row, col), J_rc) in zip(s.set.jacobian_structure, J)
        dual[col] += sign * λ[row] * J_rc
    end
    return dual
end

function MOI.get(
    model::Optimizer,
    attr::MOI.LagrangeMultiplier,
    ci::MOI.ConstraintIndex{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    λ = Float64[
        UnoSolver.uno_get_constraint_dual_solution_component(model.solver, r)
        for r in row(model, ci)
    ]
    return _dual_multiplier(model) * λ
end

function MOI.supports(
    ::Optimizer,
    ::MOI.LagrangeMultiplierStart,
    ::Type{MOI.ConstraintIndex{F,S}},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    return true
end

function MOI.get(
    model::Optimizer,
    attr::MOI.LagrangeMultiplierStart,
    ci::MOI.ConstraintIndex{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    _, cache = model.vector_nonlinear_oracle_constraints[ci.value]
    return cache.start
end

function MOI.set(
    model::Optimizer,
    attr::MOI.LagrangeMultiplierStart,
    ci::MOI.ConstraintIndex{F,S},
    start::Union{Nothing,Vector{Float64}},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    _, cache = model.vector_nonlinear_oracle_constraints[ci.value]
    cache.start = start
    return
end

### UserDefinedFunction

MOI.supports(model::Optimizer, ::MOI.UserDefinedFunction) = true

function MOI.set(model::Optimizer, attr::MOI.UserDefinedFunction, args)
    _check_no_nlp_block(model)
    MOI.Nonlinear.register_operator(
        model.model.inner,
        attr.name,
        attr.arity,
        args...,
    )
    return
end

### ListOfSupportedNonlinearOperators

function MOI.get(model::Optimizer, attr::MOI.ListOfSupportedNonlinearOperators)
    _check_no_nlp_block(model)
    return MOI.get(model.model.inner, attr)
end

### MOI.VariablePrimalStart

function MOI.supports(
    ::Optimizer,
    ::MOI.VariablePrimalStart,
    ::Type{MOI.VariableIndex},
)
    return true
end

function MOI.get(
    model::Optimizer,
    attr::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
)
    if MOI.Nonlinear._is_parameter(vi)
        throw(MOI.GetAttributeNotAllowed(attr, "Variable is a Parameter"))
    end
    MOI.throw_if_not_valid(model, vi)
    return model.variable_primal_start[column(vi)]
end

function MOI.set(
    model::Optimizer,
    attr::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
    value::Union{Real,Nothing},
)
    if MOI.Nonlinear._is_parameter(vi)
        throw(MOI.SetAttributeNotAllowed(attr, "Variable is a Parameter"))
    end
    MOI.throw_if_not_valid(model, vi)
    model.variable_primal_start[column(vi)] = value
    # No need to reset model.inner and model.solver, because this gets handled in optimize!.
    return
end

### MOI.ConstraintDualStart

_dual_start(::Optimizer, ::Nothing, ::Int = 1) = 0.0

function _dual_start(model::Optimizer, value::Real, scale::Int = 1)
    return _dual_multiplier(model) * value * scale
end

### MOI.NLPBlockDualStart

MOI.supports(::Optimizer, ::MOI.NLPBlockDualStart) = true

function MOI.set(
    model::Optimizer,
    ::MOI.NLPBlockDualStart,
    values::Union{Nothing,Vector},
)
    model.nlp_dual_start = values
    # No need to reset model.inner and model.solver, because this gets handled in optimize!.
    return
end

MOI.get(model::Optimizer, ::MOI.NLPBlockDualStart) = model.nlp_dual_start

### MOI.NLPBlock

MOI.supports(::Optimizer, ::MOI.NLPBlock) = true

# This may also be set by `optimize!` and contain the block created from
# ScalarNonlinearFunction
MOI.get(model::Optimizer, ::MOI.NLPBlock) = model.nlp_data

function MOI.set(model::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    if !MOI.is_empty(model.model.inner)
        error("Cannot mix the new and legacy nonlinear APIs")
    end
    model.nlp_data = nlp_data
    model.uses_nlp_block = !(nlp_data.evaluator isa _EmptyNLPEvaluator)
    model.inner = nothing
    model.solver = nothing
    return
end

### ObjectiveSense

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveSense,
    sense::MOI.OptimizationSense,
)
    model.sense = sense
    model.needs_new_inner = true
    model.solver = nothing
    return
end

MOI.get(model::Optimizer, ::MOI.ObjectiveSense) = model.sense

### ObjectiveFunction

function MOI.get(model::Optimizer, attr::MOI.ObjectiveFunctionType)
    if model.model.inner.objective !== nothing
        return MOI.ScalarNonlinearFunction
    end
    return MOI.get(model.model, attr)
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{
        <:Union{
            MOI.VariableIndex,
            MOI.ScalarAffineFunction{Float64},
            MOI.ScalarQuadraticFunction{Float64},
        },
    },
)
    return true
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ObjectiveFunction{F},
) where {
    F<:Union{
        MOI.VariableIndex,
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    return convert(F, MOI.get(model.model, attr))
end

function MOI.set(
    model::Optimizer,
    attr::MOI.ObjectiveFunction{F},
    func::F,
) where {
    F<:Union{
        MOI.VariableIndex,
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    # This also clears any objective of the inner nonlinear model.
    MOI.Nonlinear.set_objective(model.model, func)
    model.inner = nothing
    model.solver = nothing
    return
end

function MOI.eval_objective(model::Optimizer, x)
    if model.sense == MOI.FEASIBILITY_SENSE
        return 0.0
    elseif model.nlp_data.has_objective
        return MOI.eval_objective(model.nlp_data.evaluator, x)::Float64
    end
    return MOI.eval_objective(model.evaluator, x)
end

function MOI.eval_objective_gradient(model::Optimizer, grad, x)
    if model.sense == MOI.FEASIBILITY_SENSE
        grad .= zero(eltype(grad))
    elseif model.nlp_data.has_objective
        MOI.eval_objective_gradient(model.nlp_data.evaluator, grad, x)
    else
        MOI.eval_objective_gradient(model.evaluator, grad, x)
    end
    return
end

### _OracleNLPEvaluator

"""
    _OracleNLPEvaluator(model::Optimizer)

The inner evaluator of the `MOI.Nonlinear.EvaluatorWithQuad` of `model`: the
rows of the `MOI.VectorNonlinearOracle` constraints followed by the rows of
the `MOI.NLPBlock` evaluator.
"""
struct _OracleNLPEvaluator{E<:MOI.AbstractNLPEvaluator} <:
       MOI.AbstractNLPEvaluator
    oracles::Vector{Tuple{MOI.VectorOfVariables,_VectorNonlinearOracleCache}}
    nlp::E
    nlp_bounds::Vector{MOI.NLPBoundsPair}
    has_objective::Bool
end

function _OracleNLPEvaluator(model::Optimizer)
    return _OracleNLPEvaluator(
        model.vector_nonlinear_oracle_constraints,
        model.nlp_data.evaluator,
        model.nlp_data.constraint_bounds,
        model.nlp_data.has_objective,
    )
end

function MOI.features_available(d::_OracleNLPEvaluator)
    features = MOI.features_available(d.nlp)
    if any(s.set.eval_hessian_lagrangian === nothing for (_, s) in d.oracles)
        features = setdiff(features, [:Hess, :HessVec])
    end
    if !isempty(d.oracles)
        # The oracles do not implement the Jacobian and Hessian products.
        features = setdiff(features, [:JacVec, :HessVec])
    end
    return features
end

function MOI.initialize(d::_OracleNLPEvaluator, features::Vector{Symbol})
    return MOI.initialize(d.nlp, features)
end

MOI.eval_objective(d::_OracleNLPEvaluator, x) = MOI.eval_objective(d.nlp, x)

function MOI.eval_objective_gradient(d::_OracleNLPEvaluator, grad, x)
    return MOI.eval_objective_gradient(d.nlp, grad, x)
end

function MOI.Nonlinear._constraint_bounds(d::_OracleNLPEvaluator)
    bounds = MOI.NLPBoundsPair[]
    for (_, s) in d.oracles
        for i in 1:s.set.output_dimension
            push!(bounds, MOI.NLPBoundsPair(s.set.l[i], s.set.u[i]))
        end
    end
    return append!(bounds, d.nlp_bounds)
end

MOI.Nonlinear._has_objective(d::_OracleNLPEvaluator) = d.has_objective

function _eval_constraint(
    g::AbstractVector,
    offset::Int,
    x::AbstractVector,
    f::MOI.VectorOfVariables,
    s::_VectorNonlinearOracleCache,
)
    for i in 1:s.set.input_dimension
        s.x[i] = x[f.variables[i].value]
    end
    ret = view(g, offset .+ (1:s.set.output_dimension))
    s.eval_f_timer += @elapsed s.set.eval_f(ret, s.x)
    return offset + s.set.output_dimension
end

function MOI.eval_constraint(d::_OracleNLPEvaluator, g, x)
    offset = 0
    for (f, s) in d.oracles
        offset = _eval_constraint(g, offset, x, f, s)
    end
    MOI.eval_constraint(d.nlp, view(g, (offset+1):length(g)), x)
    return
end

function MOI.eval_constraint(model::Optimizer, g, x)
    return MOI.eval_constraint(model.evaluator, g, x)
end

function _jacobian_structure(
    ret::AbstractVector,
    row_offset::Int,
    f::MOI.VectorOfVariables,
    s::_VectorNonlinearOracleCache,
)
    for (i, j) in s.set.jacobian_structure
        push!(ret, (row_offset + i, f.variables[j].value))
    end
    return row_offset + s.set.output_dimension
end

function MOI.jacobian_structure(d::_OracleNLPEvaluator)
    J = Tuple{Int,Int}[]
    offset = 0
    for (f, s) in d.oracles
        offset = _jacobian_structure(J, offset, f, s)
    end
    if length(d.nlp_bounds) > 0
        J_nlp = MOI.jacobian_structure(d.nlp)::Vector{Tuple{Int64,Int64}}
        for (row, col) in J_nlp
            push!(J, (row + offset, col))
        end
    end
    return J
end

MOI.jacobian_structure(model::Optimizer) = MOI.jacobian_structure(model.evaluator)

function _eval_constraint_jacobian(
    values::AbstractVector,
    offset::Int,
    x::AbstractVector,
    f::MOI.VectorOfVariables,
    s::_VectorNonlinearOracleCache,
)
    for i in 1:s.set.input_dimension
        s.x[i] = x[f.variables[i].value]
    end
    nnz = length(s.set.jacobian_structure)
    s.eval_jacobian_timer +=
        @elapsed s.set.eval_jacobian(view(values, offset .+ (1:nnz)), s.x)
    return offset + nnz
end

function MOI.eval_constraint_jacobian(d::_OracleNLPEvaluator, values, x)
    offset = 0
    for (f, s) in d.oracles
        offset = _eval_constraint_jacobian(values, offset, x, f, s)
    end
    nlp_values = view(values, (offset+1):length(values))
    MOI.eval_constraint_jacobian(d.nlp, nlp_values, x)
    return
end

function MOI.eval_constraint_jacobian(model::Optimizer, values, x)
    return MOI.eval_constraint_jacobian(model.evaluator, values, x)
end

function _hessian_lagrangian_structure(
    ret::AbstractVector,
    f::MOI.VectorOfVariables,
    s::_VectorNonlinearOracleCache,
)
    for (i, j) in s.set.hessian_lagrangian_structure
        push!(ret, (f.variables[i].value, f.variables[j].value))
    end
    return
end

function MOI.hessian_lagrangian_structure(d::_OracleNLPEvaluator)
    H = Tuple{Int,Int}[]
    for (f, s) in d.oracles
        _hessian_lagrangian_structure(H, f, s)
    end
    append!(H, MOI.hessian_lagrangian_structure(d.nlp))
    return H
end

function MOI.hessian_lagrangian_structure(model::Optimizer)
    return MOI.hessian_lagrangian_structure(model.evaluator)
end

function _eval_hessian_lagrangian(
    H::AbstractVector,
    H_offset::Int,
    x::AbstractVector,
    μ::AbstractVector,
    μ_offset::Int,
    f::MOI.VectorOfVariables,
    s::_VectorNonlinearOracleCache,
)
    for i in 1:s.set.input_dimension
        s.x[i] = x[f.variables[i].value]
    end
    H_nnz = length(s.set.hessian_lagrangian_structure)
    H_view = view(H, H_offset .+ (1:H_nnz))
    μ_view = view(μ, μ_offset .+ (1:s.set.output_dimension))
    s.eval_hessian_lagrangian_timer +=
        @elapsed s.set.eval_hessian_lagrangian(H_view, s.x, μ_view)
    return H_offset + H_nnz, μ_offset + s.set.output_dimension
end

function MOI.eval_hessian_lagrangian(d::_OracleNLPEvaluator, H, x, σ, μ)
    offset, μ_offset = 0, 0
    for (f, s) in d.oracles
        offset, μ_offset =
            _eval_hessian_lagrangian(H, offset, x, μ, μ_offset, f, s)
    end
    H_nlp = view(H, (offset+1):length(H))
    μ_nlp = view(μ, (μ_offset+1):length(μ))
    MOI.eval_hessian_lagrangian(d.nlp, H_nlp, x, σ, μ_nlp)
    return
end

function MOI.eval_hessian_lagrangian(model::Optimizer, H, x, σ, μ)
    return MOI.eval_hessian_lagrangian(model.evaluator, H, x, σ, μ)
end

function MOI.eval_constraint_jacobian_product(d::_OracleNLPEvaluator, y, x, w)
    @assert isempty(d.oracles)
    return MOI.eval_constraint_jacobian_product(d.nlp, y, x, w)
end

function MOI.eval_constraint_jacobian_product(model::Optimizer, y, x, w)
    return MOI.eval_constraint_jacobian_product(model.evaluator, y, x, w)
end

function MOI.eval_constraint_jacobian_transpose_product(
    d::_OracleNLPEvaluator,
    y,
    x,
    w,
)
    @assert isempty(d.oracles)
    return MOI.eval_constraint_jacobian_transpose_product(d.nlp, y, x, w)
end

function MOI.eval_constraint_jacobian_transpose_product(model::Optimizer, y, x, w)
    return MOI.eval_constraint_jacobian_transpose_product(model.evaluator, y, x, w)
end

function MOI.eval_hessian_lagrangian_product(d::_OracleNLPEvaluator, H, x, v, σ, μ)
    @assert isempty(d.oracles)
    return MOI.eval_hessian_lagrangian_product(d.nlp, H, x, v, σ, μ)
end

function MOI.eval_hessian_lagrangian_product(model::Optimizer, H, x, v, σ, μ)
    return MOI.eval_hessian_lagrangian_product(model.evaluator, H, x, v, σ, μ)
end

### MOI.AutomaticDifferentiationBackend

MOI.supports(::Optimizer, ::MOI.AutomaticDifferentiationBackend) = true

function MOI.get(model::Optimizer, ::MOI.AutomaticDifferentiationBackend)
    return model.ad_backend
end

function MOI.set(
    model::Optimizer,
    ::MOI.AutomaticDifferentiationBackend,
    backend::MOI.Nonlinear.AbstractAutomaticDifferentiation,
)
    # Setting the backend will invalidate the model if it is different. But we
    # don't requrire == for `::MOI.Nonlinear.AutomaticDifferentiationBackend` so
    # act defensive and invalidate regardless.
    model.inner = nothing
    model.solver = nothing
    model.ad_backend = backend
    return
end

### MOI.optimize!

function _setup_model(model::Optimizer)
    vars = MOI.get(model.model.variables, MOI.ListOfVariableIndices())
    if isempty(vars)
        model.invalid_model = true
        return
    end
    # Rebuild even when the inner model is empty: a previous solve may have
    # left a stale `nlp_data` (for example, a nonlinear objective replaced by
    # a quadratic one).
    if !model.uses_nlp_block
        evaluator = MOI.Nonlinear.Evaluator(model.model.inner, model.ad_backend, vars)
        model.nlp_data = MOI.NLPBlockData(evaluator)
    end
    has_oracle = !isempty(model.vector_nonlinear_oracle_constraints)
    model.evaluator = MOI.Nonlinear.EvaluatorWithQuad(
        model.model,
        _OracleNLPEvaluator(model),
    )
    has_quadratic_constraints = any(
        isequal(MOI.Nonlinear._kFunctionTypeScalarQuadratic),
        model.model.qp.function_type,
    )
    has_nlp_constraints = !isempty(model.nlp_data.constraint_bounds) || has_oracle
    has_nlp_objective = model.nlp_data.has_objective
    features = MOI.features_available(model.evaluator)
    has_hessian = :Hess in features
    has_jacobian_operator = :JacVec in features
    has_hessian_operator = :HessVec in features
    init_feat = [:Grad]
    if has_hessian
        push!(init_feat, :Hess)
    end
    if has_hessian_operator
        push!(init_feat, :HessVec)
    end
    if has_nlp_constraints
        push!(init_feat, :Jac)
    end
    if has_jacobian_operator
        push!(init_feat, :JacVec)
    end
    MOI.initialize(model.evaluator, init_feat)

    model.hess_available = has_hessian
    model.jprod_available = has_jacobian_operator
    model.jtprod_available = has_jacobian_operator
    model.hprod_available = has_hessian_operator

    jacobian_sparsity = MOI.jacobian_structure(model)
    nnzj = length(jacobian_sparsity)
    jrows = Vector{Cint}(undef, nnzj)
    jcols = Vector{Cint}(undef, nnzj)
    for i in 1:nnzj
        jrows[i], jcols[i] = jacobian_sparsity[i]
    end
    model.jrows = jrows
    model.jcols = jcols

    hessian_sparsity = has_hessian ? MOI.hessian_lagrangian_structure(model) : Tuple{Int,Int}[]
    nnzh = length(hessian_sparsity)
    hrows = Vector{Cint}(undef, nnzh)
    hcols = Vector{Cint}(undef, nnzh)
    for i in 1:nnzh
        hrows[i], hcols[i] = hessian_sparsity[i]
    end
    model.hrows = hrows
    model.hcols = hcols

    if has_hessian && (nnzh == 0)
        model.problem_type = "LP"
    else
        if has_quadratic_constraints || has_nlp_constraints || has_nlp_objective
            model.problem_type = "NLP"
        elseif model.model.qp.objective_function_type == MOI.Nonlinear._kFunctionTypeScalarQuadratic
            model.problem_type = "QP"
        else
            @assert (model.model.qp.objective_function_type == MOI.Nonlinear._kFunctionTypeVariableIndex) || (model.model.qp.objective_function_type == MOI.Nonlinear._kFunctionTypeScalarAffine)
            model.problem_type = "LP"
        end
    end
    model.needs_new_inner = true
    return
end

function _setup_inner(model::Optimizer)::UnoSolver.Model
    if !model.needs_new_inner
        return model.inner
    end
    bounds = MOI.Nonlinear._constraint_bounds(model.evaluator)
    g_L = Float64[b.lower for b in bounds]
    g_U = Float64[b.upper for b in bounds]
    nvar = length(model.model.variables.lower)
    ncon = length(g_L)
    nnzj = length(model.jrows)
    nnzh = length(model.hrows)

    moi_objective(model, x) = MOI.eval_objective(model, x)
    moi_objective_gradient(model, g, x) = MOI.eval_objective_gradient(model, g, x)
    moi_constraints(model, c, x) = MOI.eval_constraint(model, c, x)
    moi_jacobian(model, jvals, x) = MOI.eval_constraint_jacobian(model, jvals, x)
    moi_lagrangian_hessian(model, hvals, x, multipliers, objective_multiplier) = MOI.eval_hessian_lagrangian(model, hvals, x, objective_multiplier, multipliers)
    moi_jacobian_operator(model, Jv, x, v, evaluate_at_x) = MOI.eval_constraint_jacobian_product(model, Jv, x, v)
    moi_jacobian_transposed_operator(model, Jtv, x, v, evaluate_at_x) = MOI.eval_constraint_jacobian_transpose_product(model, Jtv, x, v)
    moi_lagrangian_hessian_operator(model, Hv, x, objective_multiplier, multipliers, v, evaluate_at_x) = MOI.eval_hessian_lagrangian_product(model, Hv, x, v, objective_multiplier, multipliers)

    model.inner = UnoSolver.uno_model(
        model.problem_type,
        model.sense == MOI.MIN_SENSE,
        nvar,
        ncon,
        model.model.variables.lower,
        model.model.variables.upper,
        g_L,
        g_U,
        model.jrows,
        model.jcols,
        nnzj,
        model.hrows,
        model.hcols,
        nnzh,
        moi_objective,
        moi_constraints,
        moi_objective_gradient,
        moi_jacobian,
        model.hess_available ? moi_lagrangian_hessian : nothing,
        model.jprod_available ? moi_jacobian_operator : nothing,
        model.jtprod_available ? moi_jacobian_transposed_operator : nothing,
        model.hprod_available ? moi_lagrangian_hessian_operator : nothing,
        model,
        'L',
        1,
    )
    model.needs_new_inner = false
    return model.inner
end

function _setup_solver(model::Optimizer)::UnoSolver.Solver
    # We could reuse the C++ workspace in the future
    model.solver = UnoSolver.uno_solver()
    return model.solver
end

function MOI.optimize!(model::Optimizer)
    if model.inner === nothing
        _setup_model(model)
    end
    if model.invalid_model
        return
    end

    inner = _setup_inner(model)
    solver = _setup_solver(model)

    # The default logger is "INFO".
    UnoSolver.uno_set_solver_string_option(solver, "logger", model.silent ? "SILENT" : "INFO")

    # If provided, set the preset before any other options.
    if haskey(model.options, "preset")
        UnoSolver.uno_set_solver_preset(solver, model.options["preset"])
    end

    # Other misc options that over-ride the ones set above.
    for (name, value) in model.options
        (name == "preset") && continue
        if value isa String
            @assert UnoSolver.uno_set_solver_string_option(solver, name, value)
        elseif value isa Float64
            @assert UnoSolver.uno_set_solver_double_option(solver, name, value)
        elseif value isa Bool
            @assert UnoSolver.uno_set_solver_bool_option(solver, name, value)
        elseif value isa Integer
            @assert UnoSolver.uno_set_solver_integer_option(solver, name, Cint(value))
        else
            error(
                "Unable to add option `\"$name\"` with the value " *
                "`$value::$(typeof(value))`. The value must be a `::String`, a `::Float64`, an `::Integer`, or a `::Bool`.",
            )
        end
    end

    # Initialize the starting point, projecting variables from 0 onto their
    # bounds if VariablePrimalStart is not provided.
    for i in 1:length(model.variable_primal_start)
        x0_i = something(
            model.variable_primal_start[i],
            clamp(0.0, model.model.variables.lower[i], model.model.variables.upper[i]),
        )
        UnoSolver.uno_set_initial_primal_iterate_component(inner, i, x0_i)
    end

    for (i, start) in enumerate(model.model.qp.mult_g)
        y0_i = _dual_start(model, start, -1)
        UnoSolver.uno_set_initial_dual_iterate_component(inner, i, y0_i)
    end
    offset = length(model.model.qp.mult_g)
    if model.nlp_dual_start === nothing
        for i in offset+1:inner.ncon
            UnoSolver.uno_set_initial_dual_iterate_component(inner, i, 0.0)
        end
        # First there is VectorNonlinearOracle...
        for (_, cache) in model.vector_nonlinear_oracle_constraints
            if cache.start !== nothing
                for i in 1:cache.set.output_dimension
                    UnoSolver.uno_set_initial_dual_iterate_component(inner, offset+i, _dual_start(model, cache.start[i], -1))
                end
            end
            offset += cache.set.output_dimension
        end
        # ...then come the ScalarNonlinearFunctions
        for (key, val) in model.mult_g_nlp
            UnoSolver.uno_set_initial_dual_iterate_component(inner, offset+key.value, _dual_start(model, val, -1))
        end
    else
        for (i, start) in enumerate(model.nlp_dual_start::Vector{Float64})
            y0_i = _dual_start(model, start, -1)
            UnoSolver.uno_set_initial_dual_iterate_component(inner, offset+i, y0_i)
        end
    end

    # Clear timers
    for (_, s) in model.vector_nonlinear_oracle_constraints
        s.eval_f_timer = 0.0
        s.eval_jacobian_timer = 0.0
        s.eval_hessian_lagrangian_timer = 0.0
    end
    UnoSolver.uno_optimize(solver, inner)
    return
end

function _status_code_mapping(uno_termination_status::Cint, uno_solution_status::Cint)
    if uno_termination_status == UnoSolver.UNO_ITERATION_LIMIT
        return (MOI.ITERATION_LIMIT, MOI.UNKNOWN_RESULT_STATUS) # here we could test feasibility
    elseif uno_termination_status == UnoSolver.UNO_TIME_LIMIT
        return (MOI.TIME_LIMIT, MOI.UNKNOWN_RESULT_STATUS) # here we could test feasibility
    elseif uno_termination_status == UnoSolver.UNO_EVALUATION_ERROR
        return (MOI.INVALID_MODEL, MOI.UNKNOWN_RESULT_STATUS)
    elseif uno_termination_status == UnoSolver.UNO_ALGORITHMIC_ERROR
        return (MOI.OTHER_ERROR, MOI.UNKNOWN_RESULT_STATUS)
    elseif uno_termination_status == UnoSolver.UNO_USER_TERMINATION
        return (MOI.INTERRUPTED, MOI.UNKNOWN_RESULT_STATUS) # here we could test feasibility
    else # UNO_SUCCESS
        if uno_solution_status == UnoSolver.UNO_FEASIBLE_KKT_POINT
            return (MOI.LOCALLY_SOLVED, MOI.FEASIBLE_POINT)
        elseif uno_solution_status == UnoSolver.UNO_FEASIBLE_FJ_POINT
            return (MOI.LOCALLY_SOLVED, MOI.FEASIBLE_POINT)
        elseif uno_solution_status == UnoSolver.UNO_INFEASIBLE_STATIONARY_POINT
            return (MOI.LOCALLY_INFEASIBLE, MOI.INFEASIBLE_POINT)
        elseif uno_solution_status == UnoSolver.UNO_FEASIBLE_SMALL_STEP
            return (MOI.SLOW_PROGRESS, MOI.FEASIBLE_POINT)
        elseif uno_solution_status == UnoSolver.UNO_INFEASIBLE_SMALL_STEP
            return (MOI.SLOW_PROGRESS, MOI.INFEASIBLE_POINT)
        elseif uno_solution_status == UnoSolver.UNO_UNBOUNDED
            return (MOI.DUAL_INFEASIBLE, MOI.FEASIBLE_POINT)
        else # UNO_NOT_OPTIMAL
            return (MOI.SLOW_PROGRESS, MOI.UNKNOWN_RESULT_STATUS) # here we could test feasibility
        end
    end
end

### MOI.ResultCount

# Uno always has an iterate available.
function MOI.get(model::Optimizer, ::MOI.ResultCount)
    return (model.inner !== nothing && model.solver !== nothing) ? 1 : 0
end

### MOI.TerminationStatus

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.invalid_model
        return MOI.INVALID_MODEL
    elseif model.inner === nothing || model.solver === nothing
        return MOI.OPTIMIZE_NOT_CALLED
    end
    uno_termination_status = UnoSolver.uno_get_optimization_status(model.solver)
    uno_solution_status = UnoSolver.uno_get_solution_status(model.solver)
    termination_status, _ = _status_code_mapping(uno_termination_status, uno_solution_status)
    return termination_status
end

### MOI.RawStatusString

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    if model.invalid_model
        return "The model has no variable"
    elseif model.inner === nothing || model.solver === nothing
        return "Optimize not called"
    end
    uno_status = UnoSolver.uno_get_optimization_status(model.solver)
    return string(uno_status)
end

### MOI.PrimalStatus

function _manually_evaluated_primal_status(model::Optimizer)
    # Alexis -- revisit this code when we know how to handle the workspace
    # x, g = model.inner.x, model.inner.g
    x = Vector{Float64}(undef, model.inner.nvar)
    UnoSolver.uno_get_primal_solution(model.solver, x)
    g = Vector{Float64}(undef, model.inner.ncon)
    MOI.eval_constraint(model, g, x)

    x_L, x_U = model.model.variables.lower, model.model.variables.upper
    bounds = MOI.Nonlinear._constraint_bounds(model.evaluator)
    g_L = Float64[b.lower for b in bounds]
    g_U = Float64[b.upper for b in bounds]
    m, n = length(g_L), length(x)
    # 1e-8 is the default primal tolerance
    tol = get(model.options, "primal_tolerance", 1e-8)
    if all(x_L[i] - tol <= x[i] <= x_U[i] + tol for i in 1:n) &&
       all(g_L[i] - tol <= g[i] <= g_U[i] + tol for i in 1:m)
        return MOI.FEASIBLE_POINT
    end
    # 1e-6 is the default acceptable tolerance
    atol = get(model.options, "loose_primal_tolerance", 1e-6)
    if all(x_L[i] - atol <= x[i] <= x_U[i] + atol for i in 1:n) &&
       all(g_L[i] - atol <= g[i] <= g_U[i] + atol for i in 1:m)
        return MOI.NEARLY_FEASIBLE_POINT
    end
    return MOI.INFEASIBLE_POINT
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if !(1 <= attr.result_index <= MOI.get(model, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end
    uno_termination_status = UnoSolver.uno_get_optimization_status(model.solver)
    uno_solution_status = UnoSolver.uno_get_solution_status(model.solver)
    _, primal_status = _status_code_mapping(uno_termination_status, uno_solution_status)
    if primal_status == MOI.UNKNOWN_RESULT_STATUS
        return _manually_evaluated_primal_status(model)
    end
    return primal_status
end

### MOI.DualStatus

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    if !(1 <= attr.result_index <= MOI.get(model, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end
    uno_termination_status = UnoSolver.uno_get_optimization_status(model.solver)
    uno_solution_status = UnoSolver.uno_get_solution_status(model.solver)
    _, dual_status = _status_code_mapping(uno_termination_status, uno_solution_status)
    return dual_status
end

### MOI.SolveTimeSec

MOI.get(model::Optimizer, ::MOI.SolveTimeSec) = isnothing(model.solver) ? NaN : UnoSolver.uno_get_cpu_time(model.solver)

### MOI.BarrierIterations

MOI.get(model::Optimizer, ::MOI.BarrierIterations) = isnothing(model.solver) ? 0 : UnoSolver.uno_get_number_iterations(model.solver)

### MOI.ObjectiveValue

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    obj_val = UnoSolver.uno_get_solution_objective(model.solver)
    return obj_val
end

### MOI.VariablePrimal

function MOI.get(
    model::Optimizer,
    attr::MOI.VariablePrimal,
    vi::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, vi)
    if MOI.Nonlinear._is_parameter(vi)
        ci = MOI.ConstraintIndex{MOI.VariableIndex,MOI.Parameter{Float64}}(
            vi.value,
        )
        return MOI.get(model.model, MOI.ConstraintSet(), ci).value
    end
    return UnoSolver.uno_get_primal_solution_component(model.solver, column(vi))
end

### MOI.ConstraintPrimal

function row(
    model::Optimizer,
    ci::MOI.ConstraintIndex{F},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    return ci.value
end

function row(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction},
)
    offset = length(model.model)
    for (_, s) in model.vector_nonlinear_oracle_constraints
        offset += s.set.output_dimension
    end
    return offset + ci.value
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{<:_FUNCTIONS,<:_SETS},
)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    # Alexis -- revisit this code when we know how to handle the workspace
    # model.inner.g[row(model, ci)]
    x = Vector{Float64}(undef, model.inner.nvar)
    UnoSolver.uno_get_primal_solution(model.solver, x)
    g = Vector{Float64}(undef, model.inner.ncon)
    MOI.eval_constraint(model, g, x)
    return g[row(model, ci)]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    # Alexis -- revisit this code when we know how to handle the workspace
    # model.inner.x[ci.value]
    x = Vector{Float64}(undef, model.inner.nvar)
    UnoSolver.uno_get_primal_solution(model.solver, x)
    return x[ci.value]
end

### MOI.ConstraintDual

function _dual_multiplier(model::Optimizer)
    multiplier = (model.sense == MOI.MAX_SENSE) ? 1.0 : -1.0
    return multiplier
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{<:_FUNCTIONS,<:_SETS},
)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    λ = UnoSolver.uno_get_constraint_dual_solution_component(model.solver, row(model, ci))
    return _dual_multiplier(model) * λ
end

_reduced_cost_to_dual(::Type{S}, rc) where {S} = rc
_reduced_cost_to_dual(::Type{MOI.GreaterThan{Float64}}, rc) = max(0.0, rc)
_reduced_cost_to_dual(::Type{MOI.LessThan{Float64}}, rc) = min(0.0, rc)

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:_SETS}
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    xL_i = UnoSolver.uno_get_lower_bound_dual_solution_component(model.solver, ci.value)
    xU_i = UnoSolver.uno_get_upper_bound_dual_solution_component(model.solver, ci.value)
    return _reduced_cost_to_dual(S, _dual_multiplier(model) * (xL_i + xU_i))
end

### MOI.NLPBlockDual

function MOI.get(model::Optimizer, attr::MOI.NLPBlockDual)
    MOI.check_result_index_bounds(model, attr)
    s = _dual_multiplier(model)
    return Float64[
        s * UnoSolver.uno_get_constraint_dual_solution_component(model.solver, i)
        for i in (length(model.model)+1):model.inner.ncon
    ]
end
