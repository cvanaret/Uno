# Copyright (c) 2013: Iain Dunning, Miles Lubin, and contributors
# 2025: Adapted for UnoSolver.jl by Alexis Montoison and Charlie Vanaret
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

const _PARAMETER_OFFSET = 0x00f0000000000000

_is_parameter(x::MOI.VariableIndex) = x.value >= _PARAMETER_OFFSET
_is_parameter(term::MOI.ScalarAffineTerm) = _is_parameter(term.variable)
_is_parameter(term::MOI.ScalarQuadraticTerm) = _is_parameter(term.variable_1) || _is_parameter(term.variable_2)

mutable struct _VectorNonlinearOracleCache
    set::MOI.VectorNonlinearOracle{Float64}
    x::Vector{Float64}
    eval_f_timer::Float64
    eval_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64

    function _VectorNonlinearOracleCache(set::MOI.VectorNonlinearOracle{Float64})
        return new(set, zeros(set.input_dimension), 0.0, 0.0, 0.0)
    end
end

"""
    Optimizer()

Create a new Uno optimizer.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    model::Union{Nothing,UnoSolver.Model}
    solver::Union{Nothing,UnoSolver.Solver}
    name::String
    invalid_model::Bool
    silent::Bool
    options::Dict{String,Any}
    sense::MOI.OptimizationSense
    parameters::Dict{MOI.VariableIndex,MOI.Nonlinear.ParameterIndex}
    variables::MOI.Utilities.VariablesContainer{Float64}
    list_of_variable_indices::Vector{MOI.VariableIndex}
    variable_primal_start::Vector{Union{Nothing,Float64}}
    mult_x_L::Vector{Union{Nothing,Float64}}
    mult_x_U::Vector{Union{Nothing,Float64}}
    nlp_data::MOI.NLPBlockData
    nlp_dual_start::Union{Nothing,Vector{Float64}}
    mult_g_nlp::Dict{MOI.Nonlinear.ConstraintIndex,Float64}
    qp_data::QPBlockData{Float64}
    nlp_model::Union{Nothing,MOI.Nonlinear.Model}
    ad_backend::MOI.Nonlinear.AbstractAutomaticDifferentiation
    vector_nonlinear_oracle_constraints::Vector{Tuple{MOI.VectorOfVariables,_VectorNonlinearOracleCache}}
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
            Dict{MOI.VariableIndex,Float64}(),
            MOI.Utilities.VariablesContainer{Float64}(),
            MOI.VariableIndex[],
            Union{Nothing,Float64}[],
            Union{Nothing,Float64}[],
            Union{Nothing,Float64}[],
            MOI.NLPBlockData([], _EmptyNLPEvaluator(), false),
            nothing,
            Dict{MOI.Nonlinear.ConstraintIndex,Float64}(),
            QPBlockData{Float64}(),
            nothing,
            MOI.Nonlinear.SparseReverseMode(),
            Tuple{MOI.VectorOfVariables,_VectorNonlinearOracleCache}[],
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

function MOI.empty!(optimizer::Optimizer)
    optimizer.model = nothing
    optimizer.solver = nothing
    # SKIP: optimizer.name
    optimizer.invalid_model = false
    # SKIP: optimizer.silent
    # SKIP: optimizer.options
    optimizer.sense = MOI.FEASIBILITY_SENSE
    empty!(optimizer.parameters)
    MOI.empty!(optimizer.variables)
    empty!(optimizer.list_of_variable_indices)
    empty!(optimizer.variable_primal_start)
    empty!(optimizer.mult_x_L)
    empty!(optimizer.mult_x_U)
    optimizer.nlp_data = MOI.NLPBlockData([], _EmptyNLPEvaluator(), false)
    optimizer.nlp_dual_start = nothing
    empty!(optimizer.mult_g_nlp)
    optimizer.qp_data = QPBlockData{Float64}()
    optimizer.nlp_model = nothing
    # SKIP: optimizer.ad_backend
    empty!(optimizer.vector_nonlinear_oracle_constraints)
    optimizer.problem_type = ""
    return
end

function MOI.is_empty(optimizer::Optimizer)
    return MOI.is_empty(optimizer.variables) &&
           isempty(optimizer.variable_primal_start) &&
           isempty(optimizer.mult_x_L) &&
           isempty(optimizer.mult_x_U) &&
           optimizer.nlp_data.evaluator isa _EmptyNLPEvaluator &&
           optimizer.sense == MOI.FEASIBILITY_SENSE &&
           isempty(optimizer.vector_nonlinear_oracle_constraints)
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(optimizer::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(optimizer, src)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Uno"

function _init_nlp_model(optimizer)
    if optimizer.nlp_model === nothing
        if !(optimizer.nlp_data.evaluator isa _EmptyNLPEvaluator)
            error("Cannot mix the new and legacy nonlinear APIs")
        end
        optimizer.nlp_model = MOI.Nonlinear.Model()
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
    optimizer::Optimizer,
    set::MOI.Parameter{Float64},
)
    optimizer.model = nothing
    optimizer.solver = nothing
    _init_nlp_model(optimizer)
    p = MOI.VariableIndex(_PARAMETER_OFFSET + length(optimizer.parameters))
    push!(optimizer.list_of_variable_indices, p)
    optimizer.parameters[p] =
        MOI.Nonlinear.add_parameter(optimizer.nlp_model, set.value)
    ci = MOI.ConstraintIndex{MOI.VariableIndex,typeof(set)}(p.value)
    return p, ci
end

function MOI.is_valid(
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Parameter{Float64}},
)
    p = MOI.VariableIndex(ci.value)
    return haskey(optimizer.parameters, p)
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Parameter{Float64}},
    set::MOI.Parameter{Float64},
)
    p = optimizer.parameters[MOI.VariableIndex(ci.value)]
    optimizer.nlp_model[p] = set.value
    return
end

_replace_parameters(optimizer::Optimizer, f) = f

function _replace_parameters(optimizer::Optimizer, f::MOI.VariableIndex)
    if _is_parameter(f)
        return optimizer.parameters[f]
    end
    return f
end

function _replace_parameters(optimizer::Optimizer, f::MOI.ScalarAffineFunction)
    if any(_is_parameter, f.terms)
        g = convert(MOI.ScalarNonlinearFunction, f)
        return _replace_parameters(optimizer, g)
    end
    return f
end

function _replace_parameters(optimizer::Optimizer, f::MOI.ScalarQuadraticFunction)
    if any(_is_parameter, f.affine_terms) ||
       any(_is_parameter, f.quadratic_terms)
        g = convert(MOI.ScalarNonlinearFunction, f)
        return _replace_parameters(optimizer, g)
    end
    return f
end

function _replace_parameters(optimizer::Optimizer, f::MOI.ScalarNonlinearFunction)
    for (i, arg) in enumerate(f.args)
        f.args[i] = _replace_parameters(optimizer, arg)
    end
    return f
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:Union{MOI.VariableIndex,_FUNCTIONS}},
    ::Type{<:_SETS},
)
    return true
end

### MOI.ListOfConstraintTypesPresent

_add_scalar_nonlinear_constraints(ret, ::Nothing) = nothing

function _add_scalar_nonlinear_constraints(ret, nlp_model::MOI.Nonlinear.Model)
    for v in values(nlp_optimizer.constraints)
        F, S = MOI.ScalarNonlinearFunction, typeof(v.set)
        if !((F, S) in ret)
            push!(ret, (F, S))
        end
    end
    return
end

function MOI.get(optimizer::Optimizer, attr::MOI.ListOfConstraintTypesPresent)
    ret = MOI.get(optimizer.variables, attr)
    append!(ret, MOI.get(optimizer.qp_data, attr))
    _add_scalar_nonlinear_constraints(ret, optimizer.nlp_model)
    if !isempty(optimizer.vector_nonlinear_oracle_constraints)
        push!(ret, (MOI.VectorOfVariables, MOI.VectorNonlinearOracle{Float64}))
    end
    return ret
end

### MOI.Name

MOI.supports(::Optimizer, ::MOI.Name) = true

function MOI.set(optimizer::Optimizer, ::MOI.Name, value::String)
    optimizer.name = value
    return
end

MOI.get(optimizer::Optimizer, ::MOI.Name) = optimizer.name

### MOI.Silent

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(optimizer::Optimizer, ::MOI.Silent, value)
    optimizer.silent = value
    return
end

MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent

### MOI.TimeLimitSec

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(optimizer::Optimizer, ::MOI.TimeLimitSec, value::Real)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_limit"), Float64(value))
    return
end

function MOI.set(optimizer::Optimizer, ::MOI.TimeLimitSec, ::Nothing)
    delete!(optimizer.options, "time_limit")
    return
end

function MOI.get(optimizer::Optimizer, ::MOI.TimeLimitSec)
    return get(optimizer.options, "time_limit", nothing)
end

### MOI.RawOptimizerAttribute

MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute) = true

function MOI.set(optimizer::Optimizer, p::MOI.RawOptimizerAttribute, value)
    optimizer.options[p.name] = value
    # No need to reset optimizer.model and optimizer.solver, because this gets handled in optimize!.
    return
end

function MOI.get(optimizer::Optimizer, p::MOI.RawOptimizerAttribute)
    if !haskey(optimizer.options, p.name)
        msg = "RawOptimizerAttribute with name $(p.name) is not already set."
        throw(MOI.GetAttributeNotAllowed(p, msg))
    end
    return optimizer.options[p.name]
end

### Variables

"""
    column(x::MOI.VariableIndex)

Return the column associated with a variable.
"""
column(x::MOI.VariableIndex) = x.value

function MOI.add_variable(optimizer::Optimizer)
    push!(optimizer.variable_primal_start, nothing)
    push!(optimizer.mult_x_L, nothing)
    push!(optimizer.mult_x_U, nothing)
    optimizer.model = nothing
    optimizer.solver = nothing
    x = MOI.add_variable(optimizer.variables)
    push!(optimizer.list_of_variable_indices, x)
    return x
end

function MOI.is_valid(optimizer::Optimizer, x::MOI.VariableIndex)
    if _is_parameter(x)
        return haskey(optimizer.parameters, x)
    end
    return MOI.is_valid(optimizer.variables, x)
end

function MOI.get(optimizer::Optimizer, ::MOI.ListOfVariableIndices)
    return optimizer.list_of_variable_indices
end

function MOI.get(optimizer::Optimizer, ::MOI.NumberOfVariables)
    return length(optimizer.list_of_variable_indices)
end

function MOI.is_valid(
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    return MOI.is_valid(optimizer.variables, ci)
end

function MOI.get(
    optimizer::Optimizer,
    attr::Union{
        MOI.NumberOfConstraints{MOI.VariableIndex,<:_SETS},
        MOI.ListOfConstraintIndices{MOI.VariableIndex,<:_SETS},
    },
)
    return MOI.get(optimizer.variables, attr)
end

function MOI.get(
    optimizer::Optimizer,
    attr::Union{MOI.ConstraintFunction,MOI.ConstraintSet},
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    return MOI.get(optimizer.variables, attr, c)
end

function MOI.add_constraint(optimizer::Optimizer, x::MOI.VariableIndex, set::_SETS)
    index = MOI.add_constraint(optimizer.variables, x, set)
    optimizer.model = nothing
    optimizer.solver = nothing
    return index
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
    set::S,
) where {S<:_SETS}
    MOI.set(optimizer.variables, MOI.ConstraintSet(), ci, set)
    optimizer.model = nothing
    optimizer.solver = nothing
    return
end

function MOI.delete(
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    MOI.delete(optimizer.variables, ci)
    optimizer.model = nothing
    optimizer.solver = nothing
    return
end

### ScalarAffineFunction and ScalarQuadraticFunction constraints

function MOI.is_valid(
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex{F,<:_SETS},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    return MOI.is_valid(optimizer.qp_data, ci)
end

function MOI.add_constraint(
    optimizer::Optimizer,
    func::Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
    set::_SETS,
)
    index = MOI.add_constraint(optimizer.qp_data, func, set)
    optimizer.model = nothing
    optimizer.solver = nothing
    return index
end

function MOI.get(
    optimizer::Optimizer,
    attr::Union{MOI.NumberOfConstraints{F,S},MOI.ListOfConstraintIndices{F,S}},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
    S<:_SETS,
}
    return MOI.get(optimizer.qp_data, attr)
end

function MOI.get(
    optimizer::Optimizer,
    attr::Union{MOI.ConstraintFunction,MOI.ConstraintSet},
    c::MOI.ConstraintIndex{F,<:_SETS},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    return MOI.get(optimizer.qp_data, attr, c)
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{F,S},
    set::S,
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
    S<:_SETS,
}
    MOI.set(optimizer.qp_data, MOI.ConstraintSet(), ci, set)
    optimizer.model = nothing
    optimizer.solver = nothing
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
    optimizer::Optimizer,
    attr::MOI.ConstraintDualStart,
    c::MOI.ConstraintIndex{F,<:_SETS},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    return MOI.get(optimizer.qp_data, attr, c)
end

function MOI.set(
    optimizer::Optimizer,
    attr::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{F,<:_SETS},
    value::Union{Real,Nothing},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    MOI.throw_if_not_valid(optimizer, ci)
    MOI.set(optimizer.qp_data, attr, ci, value)
    # No need to reset optimizer.model and optimizer.solver, because this gets handled in optimize!.
    return
end

### ScalarNonlinearFunction

function MOI.is_valid(
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,<:_SETS},
)
    if optimizer.nlp_model === nothing
        return false
    end
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    return MOI.is_valid(optimizer.nlp_model, index)
end

function MOI.add_constraint(
    optimizer::Optimizer,
    f::MOI.ScalarNonlinearFunction,
    s::_SETS,
)
    _init_nlp_model(optimizer)
    if !isempty(optimizer.parameters)
        _replace_parameters(optimizer, f)
    end
    index = MOI.Nonlinear.add_constraint(optimizer.nlp_model, f, s)
    optimizer.model = nothing
    optimizer.solver = nothing
    return MOI.ConstraintIndex{typeof(f),typeof(s)}(index.value)
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ListOfConstraintIndices{F,S},
) where {F<:MOI.ScalarNonlinearFunction,S<:_SETS}
    ret = MOI.ConstraintIndex{F,S}[]
    if optimizer.nlp_model === nothing
        return ret
    end
    for (k, v) in optimizer.nlp_optimizer.constraints
        if v.set isa S
            push!(ret, MOI.ConstraintIndex{F,S}(k.value))
        end
    end
    return ret
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.NumberOfConstraints{F,S},
) where {F<:MOI.ScalarNonlinearFunction,S<:_SETS}
    if optimizer.nlp_model === nothing
        return 0
    end
    return count(v.set isa S for v in values(optimizer.nlp_optimizer.constraints))
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction},
)
    return true
end

function MOI.set(
    optimizer::Optimizer,
    attr::MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction},
    func::MOI.ScalarNonlinearFunction,
)
    _init_nlp_model(optimizer)
    if !isempty(optimizer.parameters)
        _replace_parameters(optimizer, func)
    end
    MOI.Nonlinear.set_objective(optimizer.nlp_model, func)
    optimizer.model = nothing
    optimizer.solver = nothing
    return
end

function MOI.get(
    optimizer::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,<:_SETS},
)
    MOI.throw_if_not_valid(optimizer, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    return optimizer.nlp_model[index].set
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
    set::S,
) where {S<:_SETS}
    MOI.throw_if_not_valid(optimizer, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    func = optimizer.nlp_model[index].expression
    optimizer.nlp_optimizer.constraints[index] = MOI.Nonlinear.Constraint(func, set)
    optimizer.model = nothing
    optimizer.solver = nothing
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
    optimizer::Optimizer,
    attr::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,<:_SETS},
)
    MOI.throw_if_not_valid(optimizer, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    return get(optimizer.mult_g_nlp, index, nothing)
end

function MOI.set(
    optimizer::Optimizer,
    attr::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,<:_SETS},
    value::Union{Real,Nothing},
)
    MOI.throw_if_not_valid(optimizer, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    if value === nothing
        delete!(optimizer.mult_g_nlp, index)
    else
        optimizer.mult_g_nlp[index] = convert(Float64, value)
    end
    # No need to reset optimizer.model and optimizer.solver, because this gets handled in optimize!.
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
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        MOI.VectorNonlinearOracle{Float64},
    },
)
    return 1 <= ci.value <= length(optimizer.vector_nonlinear_oracle_constraints)
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ListOfConstraintIndices{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    n = length(optimizer.vector_nonlinear_oracle_constraints)
    return MOI.ConstraintIndex{F,S}.(1:n)
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.NumberOfConstraints{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    return length(optimizer.vector_nonlinear_oracle_constraints)
end

function MOI.add_constraint(
    optimizer::Optimizer,
    f::F,
    s::S,
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    optimizer.model = nothing
    optimizer.solver = nothing
    cache = _VectorNonlinearOracleCache(s)
    push!(optimizer.vector_nonlinear_oracle_constraints, (f, cache))
    n = length(optimizer.vector_nonlinear_oracle_constraints)
    return MOI.ConstraintIndex{F,S}(n)
end

function row(
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    offset = length(optimizer.qp_data)
    for i in 1:(ci.value-1)
        _, s = optimizer.vector_nonlinear_oracle_constraints[i]
        offset += s.set.output_dimension
    end
    _, s = optimizer.vector_nonlinear_oracle_constraints[ci.value]
    return offset .+ (1:s.set.output_dimension)
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    MOI.check_result_index_bounds(optimizer, attr)
    MOI.throw_if_not_valid(optimizer, ci)
    f, _ = optimizer.vector_nonlinear_oracle_constraints[ci.value]
    return MOI.get.(optimizer, MOI.VariablePrimal(attr.result_index), f.variables)
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{F,S},
) where {F<:MOI.VectorOfVariables,S<:MOI.VectorNonlinearOracle{Float64}}
    MOI.check_result_index_bounds(optimizer, attr)
    MOI.throw_if_not_valid(optimizer, ci)
    f, s = optimizer.vector_nonlinear_oracle_constraints[ci.value]
    J = zeros(length(s.set.jacobian_structure))
    x = [MOI.get(optimizer, MOI.VariablePrimal(), xi) for xi in f.variables]
    s.set.eval_jacobian(J, x)
    λ = Float64[
        UnoSolver.uno_get_constraint_dual_solution_component(optimizer.solver, r - 1)
        for r in row(optimizer, ci)
    ]
    dual = zeros(MOI.dimension(s.set))
    sign = _dual_multiplier(optimizer)
    for ((row, col), J_rc) in zip(s.set.jacobian_structure, J)
        dual[col] += sign * λ[row] * J_rc
    end
    return dual
end

### UserDefinedFunction

MOI.supports(optimizer::Optimizer, ::MOI.UserDefinedFunction) = true

function MOI.set(optimizer::Optimizer, attr::MOI.UserDefinedFunction, args)
    _init_nlp_model(optimizer)
    MOI.Nonlinear.register_operator(
        optimizer.nlp_model,
        attr.name,
        attr.arity,
        args...,
    )
    return
end

### ListOfSupportedNonlinearOperators

function MOI.get(optimizer::Optimizer, attr::MOI.ListOfSupportedNonlinearOperators)
    _init_nlp_model(optimizer)
    return MOI.get(optimizer.nlp_model, attr)
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
    optimizer::Optimizer,
    attr::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
)
    if _is_parameter(vi)
        throw(MOI.GetAttributeNotAllowed(attr, "Variable is a Parameter"))
    end
    MOI.throw_if_not_valid(optimizer, vi)
    return optimizer.variable_primal_start[column(vi)]
end

function MOI.set(
    optimizer::Optimizer,
    attr::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
    value::Union{Real,Nothing},
)
    if _is_parameter(vi)
        throw(MOI.SetAttributeNotAllowed(attr, "Variable is a Parameter"))
    end
    MOI.throw_if_not_valid(optimizer, vi)
    optimizer.variable_primal_start[column(vi)] = value
    # No need to reset optimizer.model and optimizer.solver, because this gets handled in optimize!.
    return
end

### MOI.ConstraintDualStart

_dual_start(::Optimizer, ::Nothing, ::Int = 1) = 0.0

function _dual_start(optimizer::Optimizer, value::Real, scale::Int = 1)
    return _dual_multiplier(optimizer) * value * scale
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ConstraintDualStart,
    ::Type{MOI.ConstraintIndex{MOI.VariableIndex,S}},
) where {S<:_SETS}
    return true
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
    value::Union{Real,Nothing},
)
    MOI.throw_if_not_valid(optimizer, ci)
    optimizer.mult_x_L[ci.value] = value
    # No need to reset optimizer.model and optimizer.solver, because this gets handled in optimize!.
    return
end

function MOI.get(
    optimizer::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
)
    MOI.throw_if_not_valid(optimizer, ci)
    return optimizer.mult_x_L[ci.value]
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
    value::Union{Real,Nothing},
)
    MOI.throw_if_not_valid(optimizer, ci)
    optimizer.mult_x_U[ci.value] = value
    # No need to reset optimizer.model and optimizer.solver, because this gets handled in optimize!.
    return
end

function MOI.get(
    optimizer::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
)
    MOI.throw_if_not_valid(optimizer, ci)
    return optimizer.mult_x_U[ci.value]
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
    value::Union{Real,Nothing},
) where {S<:Union{MOI.EqualTo{Float64},MOI.Interval{Float64}}}
    MOI.throw_if_not_valid(optimizer, ci)
    if value === nothing
        optimizer.mult_x_L[ci.value] = nothing
        optimizer.mult_x_U[ci.value] = nothing
    elseif value >= 0.0
        optimizer.mult_x_L[ci.value] = value
        optimizer.mult_x_U[ci.value] = 0.0
    else
        optimizer.mult_x_L[ci.value] = 0.0
        optimizer.mult_x_U[ci.value] = value
    end
    # No need to reset optimizer.model and optimizer.solver, because this gets handled in optimize!.
    return
end

function MOI.get(
    optimizer::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:Union{MOI.EqualTo{Float64},MOI.Interval{Float64}}}
    MOI.throw_if_not_valid(optimizer, ci)
    l = optimizer.mult_x_L[ci.value]
    u = optimizer.mult_x_U[ci.value]
    return (l === u === nothing) ? nothing : (l + u)
end

### MOI.NLPBlockDualStart

MOI.supports(::Optimizer, ::MOI.NLPBlockDualStart) = true

function MOI.set(
    optimizer::Optimizer,
    ::MOI.NLPBlockDualStart,
    values::Union{Nothing,Vector},
)
    optimizer.nlp_dual_start = values
    # No need to reset optimizer.model and optimizer.solver, because this gets handled in optimize!.
    return
end

MOI.get(optimizer::Optimizer, ::MOI.NLPBlockDualStart) = optimizer.nlp_dual_start

### MOI.NLPBlock

MOI.supports(::Optimizer, ::MOI.NLPBlock) = true

# This may also be set by `optimize!` and contain the block created from
# ScalarNonlinearFunction
MOI.get(optimizer::Optimizer, ::MOI.NLPBlock) = optimizer.nlp_data

function MOI.set(optimizer::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    if optimizer.nlp_model !== nothing
        error("Cannot mix the new and legacy nonlinear APIs")
    end
    optimizer.nlp_data = nlp_data
    optimizer.model = nothing
    optimizer.solver = nothing
    return
end

### ObjectiveSense

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ObjectiveSense,
    sense::MOI.OptimizationSense,
)
    optimizer.sense = sense
    optimizer.model = nothing
    optimizer.solver = nothing
    return
end

MOI.get(optimizer::Optimizer, ::MOI.ObjectiveSense) = optimizer.sense

### ObjectiveFunction

function MOI.get(optimizer::Optimizer, attr::MOI.ObjectiveFunctionType)
    if optimizer.nlp_model !== nothing && optimizer.nlp_optimizer.objective !== nothing
        return MOI.ScalarNonlinearFunction
    end
    return MOI.get(optimizer.qp_data, attr)
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
    optimizer::Optimizer,
    attr::MOI.ObjectiveFunction{F},
) where {
    F<:Union{
        MOI.VariableIndex,
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    return convert(F, MOI.get(optimizer.qp_data, attr))
end

function MOI.set(
    optimizer::Optimizer,
    attr::MOI.ObjectiveFunction{F},
    func::F,
) where {
    F<:Union{
        MOI.VariableIndex,
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    MOI.set(optimizer.qp_data, attr, func)
    if optimizer.nlp_model !== nothing
        MOI.Nonlinear.set_objective(optimizer.nlp_model, nothing)
    end
    optimizer.model = nothing
    optimizer.solver = nothing
    return
end

function MOI.eval_objective(optimizer::Optimizer, x)
    if optimizer.sense == MOI.FEASIBILITY_SENSE
        return 0.0
    elseif optimizer.nlp_data.has_objective
        return MOI.eval_objective(optimizer.nlp_data.evaluator, x)::Float64
    end
    return MOI.eval_objective(optimizer.qp_data, x)
end

function MOI.eval_objective_gradient(optimizer::Optimizer, grad, x)
    if optimizer.sense == MOI.FEASIBILITY_SENSE
        grad .= zero(eltype(grad))
    elseif optimizer.nlp_data.has_objective
        MOI.eval_objective_gradient(optimizer.nlp_data.evaluator, grad, x)
    else
        MOI.eval_objective_gradient(optimizer.qp_data, grad, x)
    end
    return
end

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

function MOI.eval_constraint(optimizer::Optimizer, g, x)
    MOI.eval_constraint(optimizer.qp_data, g, x)
    offset = length(optimizer.qp_data)
    for (f, s) in optimizer.vector_nonlinear_oracle_constraints
        offset = _eval_constraint(g, offset, x, f, s)
    end
    g_nlp = view(g, (offset+1):length(g))
    MOI.eval_constraint(optimizer.nlp_data.evaluator, g_nlp, x)
    return
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

function MOI.jacobian_structure(optimizer::Optimizer)
    J = MOI.jacobian_structure(optimizer.qp_data)
    offset = length(optimizer.qp_data)
    for (f, s) in optimizer.vector_nonlinear_oracle_constraints
        offset = _jacobian_structure(J, offset, f, s)
    end
    if length(optimizer.nlp_data.constraint_bounds) > 0
        J_nlp = MOI.jacobian_structure(
            optimizer.nlp_data.evaluator,
        )::Vector{Tuple{Int64,Int64}}
        for (row, col) in J_nlp
            push!(J, (row + offset, col))
        end
    end
    return J
end

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

function MOI.eval_constraint_jacobian(optimizer::Optimizer, values, x)
    offset = MOI.eval_constraint_jacobian(optimizer.qp_data, values, x)
    offset -= 1  # .qp_data returns one-indexed offset
    for (f, s) in optimizer.vector_nonlinear_oracle_constraints
        offset = _eval_constraint_jacobian(values, offset, x, f, s)
    end
    nlp_values = view(values, (offset+1):length(values))
    MOI.eval_constraint_jacobian(optimizer.nlp_data.evaluator, nlp_values, x)
    return
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

function MOI.hessian_lagrangian_structure(optimizer::Optimizer)
    H = MOI.hessian_lagrangian_structure(optimizer.qp_data)
    for (f, s) in optimizer.vector_nonlinear_oracle_constraints
        _hessian_lagrangian_structure(H, f, s)
    end
    append!(H, MOI.hessian_lagrangian_structure(optimizer.nlp_data.evaluator))
    return H
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

function MOI.eval_hessian_lagrangian(optimizer::Optimizer, H, x, σ, μ)
    offset = MOI.eval_hessian_lagrangian(optimizer.qp_data, H, x, σ, μ)
    offset -= 1  # optimizer.qp_data returns one-indexed offset
    μ_offset = length(optimizer.qp_data)
    for (f, s) in optimizer.vector_nonlinear_oracle_constraints
        offset, μ_offset =
            _eval_hessian_lagrangian(H, offset, x, μ, μ_offset, f, s)
    end
    H_nlp = view(H, (offset+1):length(H))
    μ_nlp = view(μ, (μ_offset+1):length(μ))
    MOI.eval_hessian_lagrangian(optimizer.nlp_data.evaluator, H_nlp, x, σ, μ_nlp)
    return
end

function MOI.eval_constraint_jacobian_product(optimizer::Optimizer, y, x, w)
    fill!(y, 0.0)
    qp_offset = length(optimizer.qp_data)
    y_nlp = view(y, (qp_offset+1):length(y))
    MOI.eval_constraint_jacobian_product(optimizer.nlp_data.evaluator, y_nlp, x, w)
    MOI.eval_constraint_jacobian_product(optimizer.qp_data, y, x, w)
    return
end

function MOI.eval_constraint_jacobian_transpose_product(optimizer::Optimizer, y, x, w)
    fill!(y, 0.0)
    qp_offset = length(optimizer.qp_data)
    w_nlp = view(w, (qp_offset+1):length(w))
    MOI.eval_constraint_jacobian_transpose_product(optimizer.nlp_data.evaluator, y, x, w_nlp)
    MOI.eval_constraint_jacobian_transpose_product(optimizer.qp_data, y, x, w)
    return
end

function MOI.eval_hessian_lagrangian_product(optimizer::Optimizer, H, x, v, σ, μ)
    fill!(H, 0.0)
    qp_offset = length(optimizer.qp_data)
    μ_nlp = view(μ, (qp_offset+1):length(μ))
    MOI.eval_hessian_lagrangian_product(optimizer.nlp_data.evaluator, H, x, v, σ, μ_nlp)
    MOI.eval_hessian_lagrangian_product(optimizer.qp_data, H, x, v, σ, μ)
    return
end

### MOI.AutomaticDifferentiationBackend

MOI.supports(::Optimizer, ::MOI.AutomaticDifferentiationBackend) = true

function MOI.get(optimizer::Optimizer, ::MOI.AutomaticDifferentiationBackend)
    return optimizer.ad_backend
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.AutomaticDifferentiationBackend,
    backend::MOI.Nonlinear.AbstractAutomaticDifferentiation,
)
    # Setting the backend will invalidate the model if it is different. But we
    # don't requrire == for `::MOI.Nonlinear.AutomaticDifferentiationBackend` so
    # act defensive and invalidate regardless.
    optimizer.model = nothing
    optimizer.solver = nothing
    optimizer.ad_backend = backend
    return
end

### MOI.optimize!

function _setup_model(optimizer::Optimizer)
    vars = MOI.get(optimizer.variables, MOI.ListOfVariableIndices())
    if isempty(vars)
        optimizer.invalid_model = true
        return
    end
    if optimizer.nlp_model !== nothing
        evaluator = MOI.Nonlinear.Evaluator(optimizer.nlp_model, optimizer.ad_backend, vars)
        optimizer.nlp_data = MOI.NLPBlockData(evaluator)
    end
    has_oracle = !isempty(optimizer.vector_nonlinear_oracle_constraints)
    has_quadratic_constraints =
        any(isequal(_kFunctionTypeScalarQuadratic), optimizer.qp_data.function_type)
    has_nlp_constraints = !isempty(optimizer.nlp_data.constraint_bounds) || has_oracle
    has_nlp_objective = optimizer.nlp_data.has_objective
    has_hessian = :Hess in MOI.features_available(optimizer.nlp_data.evaluator)
    has_jacobian_operator = :JacVec in MOI.features_available(optimizer.nlp_data.evaluator)
    has_hessian_operator = :HessVec in MOI.features_available(optimizer.nlp_data.evaluator)
    for (_, s) in optimizer.vector_nonlinear_oracle_constraints
        if s.set.eval_hessian_lagrangian === nothing
            has_hessian = false
            break
        end
    end
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
    MOI.initialize(optimizer.nlp_data.evaluator, init_feat)

    jacobian_sparsity = MOI.jacobian_structure(optimizer)
    nnzj = length(jacobian_sparsity)
    jrows = Vector{Cint}(undef, nnzj)
    jcols = Vector{Cint}(undef, nnzj)
    for i in 1:nnzj
        jrows[i], jcols[i] = jacobian_sparsity[i]
    end

    hessian_sparsity = has_hessian ? MOI.hessian_lagrangian_structure(optimizer) : Tuple{Int,Int}[]
    nnzh = length(hessian_sparsity)
    hrows = Vector{Cint}(undef, nnzh)
    hcols = Vector{Cint}(undef, nnzh)
    for i in 1:nnzh
        hrows[i], hcols[i] = hessian_sparsity[i]
    end
    if has_quadratic_constraints || has_nlp_constraints || has_nlp_objective
        optimizer.problem_type = "NLP"
    elseif optimizer.qp_data.objective_function_type == _kFunctionTypeScalarQuadratic
        optimizer.problem_type = "QP"
    else
        @assert (optimizer.qp_data.objective_function_type == _kFunctionTypeVariableIndex) || (optimizer.qp_data.objective_function_type == _kFunctionTypeScalarAffine)
        optimizer.problem_type = "LP"
    end
    if isempty(hessian_sparsity)
        optimizer.problem_type = "LP"
    end

    moi_objective(model, x) = MOI.eval_objective(model, x)
    moi_objective_gradient(model, g, x) = MOI.eval_objective_gradient(model, g, x)
    moi_constraints(model, c, x) = MOI.eval_constraint(model, c, x)
    moi_jacobian(model, jvals, x) = MOI.eval_constraint_jacobian(model, jvals, x)
    moi_lagrangian_hessian(model, hvals, x, multipliers, objective_multiplier) = MOI.eval_hessian_lagrangian(model, hvals, x, objective_multiplier, multipliers)
    moi_jacobian_operator(model, Jv, x, v, evaluate_at_x) = MOI.eval_constraint_jacobian_product(model, Jv, x, v)
    moi_jacobian_transposed_operator(model, Jtv, x, v, evaluate_at_x) = MOI.eval_constraint_jacobian_transpose_product(model, Jtv, x, v)
    moi_lagrangian_hessian_operator(model, Hv, x, objective_multiplier, multipliers, v, evaluate_at_x) = MOI.eval_hessian_lagrangian_product(model, Hv, x, v, objective_multiplier, multipliers)

    g_L, g_U = copy(optimizer.qp_data.g_L), copy(optimizer.qp_data.g_U)
    for (_, s) in optimizer.vector_nonlinear_oracle_constraints
        append!(g_L, s.set.l)
        append!(g_U, s.set.u)
    end
    for bound in optimizer.nlp_data.constraint_bounds
        push!(g_L, bound.lower)
        push!(g_U, bound.upper)
    end
    nvar = length(vars)
    ncon = length(g_L)
    optimizer.model = UnoSolver.uno_model(
        optimizer.problem_type,
        optimizer.sense == MOI.MIN_SENSE,
        nvar,
        ncon,
        optimizer.variables.lower,
        optimizer.variables.upper,
        g_L,
        g_U,
        jrows,
        jcols,
        nnzj,
        hrows,
        hcols,
        nnzh,
        moi_objective,
        moi_constraints,
        moi_objective_gradient,
        moi_jacobian,
        moi_lagrangian_hessian,
        (has_jacobian_operator && !has_oracle) ? moi_jacobian_operator : nothing,
        (has_jacobian_operator && !has_oracle) ? moi_jacobian_transposed_operator : nothing,
        (has_hessian_operator && !has_oracle) ? moi_lagrangian_hessian_operator : nothing,
        optimizer,
        'L',
        1.0,
    )
    optimizer.solver = UnoSolver.uno_solver()
    return
end

function copy_parameters(optimizer::Optimizer)
    if optimizer.nlp_model === nothing
        return
    end
    empty!(optimizer.qp_data.parameters)
    for (p, index) in optimizer.parameters
        optimizer.qp_data.parameters[p.value] = optimizer.nlp_model[index]
    end
    return
end

function MOI.optimize!(optimizer::Optimizer)
    if optimizer.model === nothing || optimizer.solver == nothing
        _setup_model(optimizer)
    end
    if optimizer.invalid_model
        return
    end
    copy_parameters(optimizer)
    model = optimizer.model::UnoSolver.Model
    solver = optimizer.solver::UnoSolver.Solver

    # The default logger is "INFO".
    UnoSolver.uno_set_solver_string_option(solver, "logger", optimizer.silent ? "SILENT" : "INFO")

    # Other misc options that over-ride the ones set above.
    for (name, value) in optimizer.options
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
    for i in 1:length(optimizer.variable_primal_start)
        x0_i = if optimizer.variable_primal_start[i] !== nothing
            optimizer.variable_primal_start[i]
        else
            clamp(0.0, optimizer.variables.lower[i], optimizer.variables.upper[i])
        end
        UnoSolver.uno_set_initial_primal_iterate_component(model, i-1, x0_i)
    end

    for (i, start) in enumerate(optimizer.qp_data.mult_g)
        y0_i = _dual_start(optimizer, start, -1)
        UnoSolver.uno_set_initial_dual_iterate_component(model, i-1, y0_i)
    end
    offset = length(optimizer.qp_data.mult_g)
    if optimizer.nlp_dual_start === nothing
        for i in offset+1:model.ncon
            UnoSolver.uno_set_initial_dual_iterate_component(model, i-1, 0.0)
        end
        for (key, val) in optimizer.mult_g_nlp
            UnoSolver.uno_set_initial_dual_iterate_component(model, offset+key.value-1, val)
        end
    else
        for (i, start) in enumerate(optimizer.nlp_dual_start::Vector{Float64})
            y0_i = _dual_start(optimizer, start, -1)
            UnoSolver.uno_set_initial_dual_iterate_component(model, offset+i-1, y0_i)
        end
    end

    # Clear timers
    for (_, s) in optimizer.vector_nonlinear_oracle_constraints
        s.eval_f_timer = 0.0
        s.eval_jacobian_timer = 0.0
        s.eval_hessian_lagrangian_timer = 0.0
    end
    UnoSolver.uno_optimize(solver, model)
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
        else # UNO_UNBOUNDED
            return (MOI.DUAL_INFEASIBLE, MOI.FEASIBLE_POINT)
        end
    end
end

### MOI.ResultCount

# Uno always has an iterate available.
function MOI.get(optimizer::Optimizer, ::MOI.ResultCount)
    return (optimizer.model !== nothing && optimizer.solver !== nothing) ? 1 : 0
end

### MOI.TerminationStatus

function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    if optimizer.invalid_model
        return MOI.INVALID_MODEL
    elseif optimizer.model === nothing || optimizer.solver === nothing
        return MOI.OPTIMIZE_NOT_CALLED
    end
    uno_termination_status = UnoSolver.uno_get_optimization_status(optimizer.solver)
    uno_solution_status = UnoSolver.uno_get_solution_status(optimizer.solver)
    termination_status, _ = _status_code_mapping(uno_termination_status, uno_solution_status)
    return termination_status
end

### MOI.RawStatusString

function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
    if optimizer.invalid_model
        return "The model has no variable"
    elseif optimizer.model === nothing || optimizer.solver === nothing
        return "Optimize not called"
    end
    uno_status = UnoSolver.uno_get_optimization_status(optimizer.solver)
    return string(uno_status)
end

### MOI.PrimalStatus

function _manually_evaluated_primal_status(optimizer::Optimizer)
    # Alexis -- revisit this code when we know how to handle the workspace
    # x, g = optimizer.model.x, optimizer.model.g
    x = Vector{Float64}(undef, optimizer.model.nvar)
    UnoSolver.uno_get_primal_solution(optimizer.solver, x)
    g = Vector{Float64}(undef, optimizer.model.ncon)
    MOI.eval_constraint(optimizer, g, x)

    x_L, x_U = optimizer.variables.lower, optimizer.variables.upper
    g_L, g_U = copy(optimizer.qp_data.g_L), copy(optimizer.qp_data.g_U)
    # Assuming constraints are guaranteed to be in the order:
    # [qp_cons, nlp_cons, oracle]
    for bound in optimizer.nlp_data.constraint_bounds
        push!(g_L, bound.lower)
        push!(g_U, bound.upper)
    end
    for (_, cache) in optimizer.vector_nonlinear_oracle_constraints
        append!(g_L, cache.set.l)
        append!(g_U, cache.set.u)
    end
    m, n = length(g_L), length(x)
    # 1e-8 is the default primal tolerance
    tol = get(optimizer.options, "primal_tolerance", 1e-8)
    if all(x_L[i] - tol <= x[i] <= x_U[i] + tol for i in 1:n) &&
       all(g_L[i] - tol <= g[i] <= g_U[i] + tol for i in 1:m)
        return MOI.FEASIBLE_POINT
    end
    # 1e-6 is the default acceptable tolerance
    atol = get(optimizer.options, "loose_primal_tolerance", 1e-6)
    if all(x_L[i] - atol <= x[i] <= x_U[i] + atol for i in 1:n) &&
       all(g_L[i] - atol <= g[i] <= g_U[i] + atol for i in 1:m)
        return MOI.NEARLY_FEASIBLE_POINT
    end
    return MOI.INFEASIBLE_POINT
end

function MOI.get(optimizer::Optimizer, attr::MOI.PrimalStatus)
    if !(1 <= attr.result_index <= MOI.get(optimizer, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end
    uno_termination_status = UnoSolver.uno_get_optimization_status(optimizer.solver)
    uno_solution_status = UnoSolver.uno_get_solution_status(optimizer.solver)
    _, primal_status = _status_code_mapping(uno_termination_status, uno_solution_status)
    if primal_status == MOI.UNKNOWN_RESULT_STATUS
        return _manually_evaluated_primal_status(optimizer)
    end
    return primal_status
end

### MOI.DualStatus

function MOI.get(optimizer::Optimizer, attr::MOI.DualStatus)
    if !(1 <= attr.result_index <= MOI.get(optimizer, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end
    uno_termination_status = UnoSolver.uno_get_optimization_status(optimizer.solver)
    uno_solution_status = UnoSolver.uno_get_solution_status(optimizer.solver)
    _, dual_status = _status_code_mapping(uno_termination_status, uno_solution_status)
    return dual_status
end

### MOI.SolveTimeSec

MOI.get(optimizer::Optimizer, ::MOI.SolveTimeSec) = isnothing(optimizer.solver) ? NaN : UnoSolver.uno_get_cpu_time(optimizer.solver)

### MOI.BarrierIterations

MOI.get(optimizer::Optimizer, ::MOI.BarrierIterations) = isnothing(optimizer.solver) ? 0 : UnoSolver.uno_get_number_iterations(optimizer.solver)

### MOI.ObjectiveValue

function MOI.get(optimizer::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    obj_val = UnoSolver.uno_get_solution_objective(optimizer.solver)
    return obj_val
end

### MOI.VariablePrimal

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.VariablePrimal,
    vi::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(optimizer, attr)
    MOI.throw_if_not_valid(optimizer, vi)
    if _is_parameter(vi)
        p = optimizer.parameters[vi]
        return optimizer.nlp_model[p]
    end
    return UnoSolver.uno_get_primal_solution_component(optimizer.solver, column(vi)-1)
end

### MOI.ConstraintPrimal

function row(
    optimizer::Optimizer,
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
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction},
)
    offset = length(optimizer.qp_data)
    for (_, s) in optimizer.vector_nonlinear_oracle_constraints
        offset += s.set.output_dimension
    end
    return offset + ci.value
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{<:_FUNCTIONS,<:_SETS},
)
    MOI.check_result_index_bounds(optimizer, attr)
    MOI.throw_if_not_valid(optimizer, ci)
    # Alexis -- revisit this code when we know how to handle the workspace
    # optimizer.model.g[row(optimizer, ci)]
    x = Vector{Float64}(undef, optimizer.model.nvar)
    UnoSolver.uno_get_primal_solution(optimizer.solver, x)
    g = Vector{Float64}(undef, optimizer.model.ncon)
    MOI.eval_constraint(optimizer, g, x)
    return g[row(optimizer, ci)]
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    MOI.check_result_index_bounds(optimizer, attr)
    MOI.throw_if_not_valid(optimizer, ci)
    # Alexis -- revisit this code when we know how to handle the workspace
    # optimizer.model.x[ci.value]
    x = Vector{Float64}(undef, optimizer.model.nvar)
    UnoSolver.uno_get_primal_solution(optimizer.solver, x)
    return x[ci.value]
end

### MOI.ConstraintDual

function _dual_multiplier(optimizer::Optimizer)
    return xor(optimizer.problem_type == "LP", optimizer.sense == MOI.MAX_SENSE) ? 1.0 : -1.0
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{<:_FUNCTIONS,<:_SETS},
)
    MOI.check_result_index_bounds(optimizer, attr)
    MOI.throw_if_not_valid(optimizer, ci)
    λ = UnoSolver.uno_get_constraint_dual_solution_component(optimizer.solver, row(optimizer, ci)-1)
    return _dual_multiplier(optimizer) * λ
end

_reduced_cost_to_dual(::Type{S}, rc) where {S} = rc
_reduced_cost_to_dual(::Type{MOI.GreaterThan{Float64}}, rc) = max(0.0, rc)
_reduced_cost_to_dual(::Type{MOI.LessThan{Float64}}, rc) = min(0.0, rc)

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:_SETS}
    MOI.check_result_index_bounds(optimizer, attr)
    MOI.throw_if_not_valid(optimizer, ci)
    xL_i = UnoSolver.uno_get_lower_bound_dual_solution_component(optimizer.solver, ci.value-1)
    xU_i = UnoSolver.uno_get_upper_bound_dual_solution_component(optimizer.solver, ci.value-1)
    return _reduced_cost_to_dual(S, _dual_multiplier(optimizer) * (xL_i + xU_i))
end

### MOI.NLPBlockDual

function MOI.get(optimizer::Optimizer, attr::MOI.NLPBlockDual)
    MOI.check_result_index_bounds(optimizer, attr)
    s = _dual_multiplier(optimizer)
    return Float64[
        s * UnoSolver.uno_get_constraint_dual_solution_component(optimizer.solver, i)
        for i in length(optimizer.qp_data):(optimizer.model.ncon-1)
    ]
end
