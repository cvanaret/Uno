# Copyright (c) 2013: Iain Dunning, Miles Lubin, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# !!! warning
#
#     The contents of this file are experimental.
#
#     Until this message is removed, breaking changes to the functions and
#     types, including their deletion, may be introduced in any minor or patch
#     release of Ipopt.

@enum(
    _FunctionType,
    _kFunctionTypeVariableIndex,
    _kFunctionTypeScalarAffine,
    _kFunctionTypeScalarQuadratic,
)

function _function_type_to_func(::Type{T}, k::_FunctionType) where {T}
    if k == _kFunctionTypeVariableIndex
        return MOI.VariableIndex
    elseif k == _kFunctionTypeScalarAffine
        return MOI.ScalarAffineFunction{T}
    else
        @assert k == _kFunctionTypeScalarQuadratic
        return MOI.ScalarQuadraticFunction{T}
    end
end

_function_info(::MOI.VariableIndex) = _kFunctionTypeVariableIndex
_function_info(::MOI.ScalarAffineFunction) = _kFunctionTypeScalarAffine
_function_info(::MOI.ScalarQuadraticFunction) = _kFunctionTypeScalarQuadratic

@enum(
    _BoundType,
    _kBoundTypeLessThan,
    _kBoundTypeGreaterThan,
    _kBoundTypeEqualTo,
    _kBoundTypeInterval,
)

_set_info(s::MOI.LessThan) = _kBoundTypeLessThan, -Inf, s.upper
_set_info(s::MOI.GreaterThan) = _kBoundTypeGreaterThan, s.lower, Inf
_set_info(s::MOI.EqualTo) = _kBoundTypeEqualTo, s.value, s.value
_set_info(s::MOI.Interval) = _kBoundTypeInterval, s.lower, s.upper

function _bound_type_to_set(::Type{T}, k::_BoundType) where {T}
    if k == _kBoundTypeEqualTo
        return MOI.EqualTo{T}
    elseif k == _kBoundTypeLessThan
        return MOI.LessThan{T}
    elseif k == _kBoundTypeGreaterThan
        return MOI.GreaterThan{T}
    else
        @assert k == _kBoundTypeInterval
        return MOI.Interval{T}
    end
end

mutable struct QPBlockData{T}
    objective::Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}}
    objective_function_type::_FunctionType
    constraints::Vector{
        Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}},
    }
    g_L::Vector{T}
    g_U::Vector{T}
    mult_g::Vector{Union{Nothing,T}}
    function_type::Vector{_FunctionType}
    bound_type::Vector{_BoundType}
    parameters::Dict{Int64,T}

    function QPBlockData{T}() where {T}
        return new(
            zero(MOI.ScalarQuadraticFunction{T}),
            _kFunctionTypeScalarAffine,
            Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}}[],
            T[],
            T[],
            Union{Nothing,T}[],
            _FunctionType[],
            _BoundType[],
            Dict{Int64,T}(),
        )
    end
end

function _value(variable::MOI.VariableIndex, x::Vector, p::Dict)
    if _is_parameter(variable)
        return p[variable.value]
    else
        return x[variable.value]
    end
end

function eval_function(
    f::MOI.ScalarQuadraticFunction{T},
    x::Vector{T},
    p::Dict{Int64,T},
)::T where {T}
    y = f.constant
    for term in f.affine_terms
        y += term.coefficient * _value(term.variable, x, p)
    end
    for term in f.quadratic_terms
        v1 = _value(term.variable_1, x, p)
        v2 = _value(term.variable_2, x, p)
        if term.variable_1 == term.variable_2
            y += term.coefficient * v1 * v2 / 2
        else
            y += term.coefficient * v1 * v2
        end
    end
    return y
end

function eval_function(
    f::MOI.ScalarAffineFunction{T},
    x::Vector{T},
    p::Dict{Int64,T},
)::T where {T}
    y = f.constant
    for term in f.terms
        y += term.coefficient * _value(term.variable, x, p)
    end
    return y
end

function eval_dense_gradient(
    ∇f::Vector{T},
    f::MOI.ScalarQuadraticFunction{T},
    x::Vector{T},
    p::Dict{Int64,T},
)::Nothing where {T}
    for term in f.affine_terms
        if !_is_parameter(term.variable)
            ∇f[term.variable.value] += term.coefficient
        end
    end
    for term in f.quadratic_terms
        if !_is_parameter(term.variable_1)
            v = _value(term.variable_2, x, p)
            ∇f[term.variable_1.value] += term.coefficient * v
        end
        if term.variable_1 != term.variable_2 && !_is_parameter(term.variable_2)
            v = _value(term.variable_1, x, p)
            ∇f[term.variable_2.value] += term.coefficient * v
        end
    end
    return
end

function eval_dense_gradient(
    ∇f::Vector{T},
    f::MOI.ScalarAffineFunction{T},
    x::Vector{T},
    p::Dict{Int64,T},
)::Nothing where {T}
    for term in f.terms
        if !_is_parameter(term.variable)
            ∇f[term.variable.value] += term.coefficient
        end
    end
    return
end

function append_sparse_gradient_structure!(
    f::MOI.ScalarQuadraticFunction,
    J,
    row,
)
    for term in f.affine_terms
        if !_is_parameter(term.variable)
            push!(J, (row, term.variable.value))
        end
    end
    for term in f.quadratic_terms
        if !_is_parameter(term.variable_1)
            push!(J, (row, term.variable_1.value))
        end
        if term.variable_1 != term.variable_2 && !_is_parameter(term.variable_2)
            push!(J, (row, term.variable_2.value))
        end
    end
    return
end

function append_sparse_gradient_structure!(f::MOI.ScalarAffineFunction, J, row)
    for term in f.terms
        if !_is_parameter(term.variable)
            push!(J, (row, term.variable.value))
        end
    end
    return
end

function eval_sparse_gradient(
    ∇f::AbstractVector{T},
    f::MOI.ScalarQuadraticFunction{T},
    x::Vector{T},
    p::Dict{Int64,T},
)::Int where {T}
    i = 0
    for term in f.affine_terms
        if !_is_parameter(term.variable)
            i += 1
            ∇f[i] = term.coefficient
        end
    end
    for term in f.quadratic_terms
        if !_is_parameter(term.variable_1)
            v = _value(term.variable_2, x, p)
            i += 1
            ∇f[i] = term.coefficient * v
        end
        if term.variable_1 != term.variable_2 && !_is_parameter(term.variable_2)
            v = _value(term.variable_1, x, p)
            i += 1
            ∇f[i] = term.coefficient * v
        end
    end
    return i
end

function eval_sparse_gradient(
    ∇f::AbstractVector{T},
    f::MOI.ScalarAffineFunction{T},
    x::Vector{T},
    p::Dict{Int64,T},
)::Int where {T}
    i = 0
    for term in f.terms
        if !_is_parameter(term.variable)
            i += 1
            ∇f[i] = term.coefficient
        end
    end
    return i
end

function append_sparse_hessian_structure!(f::MOI.ScalarQuadraticFunction, H)
    for term in f.quadratic_terms
        if _is_parameter(term.variable_1) || _is_parameter(term.variable_2)
            continue
        end
        push!(H, (term.variable_1.value, term.variable_2.value))
    end
    return
end

append_sparse_hessian_structure!(::MOI.ScalarAffineFunction, H) = nothing

function eval_sparse_hessian(
    ∇²f::AbstractVector{T},
    f::MOI.ScalarQuadraticFunction{T},
    σ::T,
)::Int where {T}
    i = 0
    for term in f.quadratic_terms
        if _is_parameter(term.variable_1) || _is_parameter(term.variable_2)
            continue
        end
        i += 1
        ∇²f[i] = term.coefficient * σ
    end
    return i
end

function eval_sparse_hessian(
    ∇²f::AbstractVector{T},
    f::MOI.ScalarAffineFunction{T},
    σ::T,
)::Int where {T}
    return 0
end

Base.length(block::QPBlockData) = length(block.bound_type)

function MOI.set(
    block::QPBlockData{T},
    ::MOI.ObjectiveFunction{F},
    f::F,
) where {T,F<:Union{MOI.VariableIndex,MOI.ScalarAffineFunction{T}}}
    block.objective = convert(MOI.ScalarAffineFunction{T}, f)
    block.objective_function_type = _function_info(f)
    return
end

function MOI.set(
    block::QPBlockData{T},
    ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}},
    f::MOI.ScalarQuadraticFunction{T},
) where {T}
    block.objective = f
    block.objective_function_type = _function_info(f)
    return
end

function MOI.get(block::QPBlockData{T}, ::MOI.ObjectiveFunctionType) where {T}
    return _function_type_to_func(T, block.objective_function_type)
end

function MOI.get(block::QPBlockData{T}, ::MOI.ObjectiveFunction{F}) where {T,F}
    return convert(F, block.objective)
end

function MOI.get(
    block::QPBlockData{T},
    ::MOI.ListOfConstraintTypesPresent,
) where {T}
    constraints = Set{Tuple{Type,Type}}()
    for i in 1:length(block)
        F = _function_type_to_func(T, block.function_type[i])
        S = _bound_type_to_set(T, block.bound_type[i])
        push!(constraints, (F, S))
    end
    return collect(constraints)
end

function MOI.is_valid(
    block::QPBlockData{T},
    ci::MOI.ConstraintIndex{F,S},
) where {
    T,
    F<:Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}},
    S<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T},MOI.Interval{T}},
}
    return 1 <= ci.value <= length(block)
end

function MOI.get(
    block::QPBlockData{T},
    ::MOI.ListOfConstraintIndices{F,S},
) where {
    T,
    F<:Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}},
    S<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T},MOI.Interval{T}},
}
    ret = MOI.ConstraintIndex{F,S}[]
    for i in 1:length(block)
        if _bound_type_to_set(T, block.bound_type[i]) != S
            continue
        elseif _function_type_to_func(T, block.function_type[i]) != F
            continue
        end
        push!(ret, MOI.ConstraintIndex{F,S}(i))
    end
    return ret
end

function MOI.get(
    block::QPBlockData{T},
    ::MOI.NumberOfConstraints{F,S},
) where {
    T,
    F<:Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}},
    S<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T},MOI.Interval{T}},
}
    return length(MOI.get(block, MOI.ListOfConstraintIndices{F,S}()))
end

function MOI.add_constraint(
    block::QPBlockData{T},
    f::Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}},
    s::Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T},MOI.Interval{T}},
) where {T}
    push!(block.constraints, f)
    bound_type, l, u = _set_info(s)
    push!(block.g_L, l)
    push!(block.g_U, u)
    push!(block.mult_g, nothing)
    push!(block.bound_type, bound_type)
    push!(block.function_type, _function_info(f))
    return MOI.ConstraintIndex{typeof(f),typeof(s)}(length(block.bound_type))
end

function MOI.get(
    block::QPBlockData{T},
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{F,S},
) where {T,F,S}
    return convert(F, block.constraints[c.value])
end

function MOI.get(
    block::QPBlockData{T},
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{F,S},
) where {T,F,S}
    row = c.value
    if block.bound_type[row] == _kBoundTypeEqualTo
        return MOI.EqualTo(block.g_L[row])
    elseif block.bound_type[row] == _kBoundTypeLessThan
        return MOI.LessThan(block.g_U[row])
    elseif block.bound_type[row] == _kBoundTypeGreaterThan
        return MOI.GreaterThan(block.g_L[row])
    else
        @assert block.bound_type[row] == _kBoundTypeInterval
        return MOI.Interval(block.g_L[row], block.g_U[row])
    end
end

function MOI.set(
    block::QPBlockData{T},
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{F,MOI.LessThan{T}},
    set::MOI.LessThan{T},
) where {T,F}
    row = c.value
    block.g_U[row] = set.upper
    return
end

function MOI.set(
    block::QPBlockData{T},
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{F,MOI.GreaterThan{T}},
    set::MOI.GreaterThan{T},
) where {T,F}
    row = c.value
    block.g_L[row] = set.lower
    return
end

function MOI.set(
    block::QPBlockData{T},
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{F,MOI.EqualTo{T}},
    set::MOI.EqualTo{T},
) where {T,F}
    row = c.value
    block.g_L[row] = set.value
    block.g_U[row] = set.value
    return
end

function MOI.set(
    block::QPBlockData{T},
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{F,MOI.Interval{T}},
    set::MOI.Interval{T},
) where {T,F}
    row = c.value
    block.g_L[row] = set.lower
    block.g_U[row] = set.upper
    return
end

function MOI.get(
    block::QPBlockData{T},
    ::MOI.ConstraintDualStart,
    c::MOI.ConstraintIndex{F,S},
) where {T,F,S}
    return block.mult_g[c.value]
end

function MOI.set(
    block::QPBlockData{T},
    ::MOI.ConstraintDualStart,
    c::MOI.ConstraintIndex{F,S},
    value,
) where {T,F,S}
    block.mult_g[c.value] = value
    return
end

function MOI.eval_objective(
    block::QPBlockData{T},
    x::AbstractVector{T},
) where {T}
    return eval_function(block.objective, x, block.parameters)
end

function MOI.eval_objective_gradient(
    block::QPBlockData{T},
    ∇f::AbstractVector{T},
    x::AbstractVector{T},
) where {T}
    ∇f .= zero(T)
    eval_dense_gradient(∇f, block.objective, x, block.parameters)
    return
end

function MOI.eval_constraint(
    block::QPBlockData{T},
    g::AbstractVector{T},
    x::AbstractVector{T},
) where {T}
    for (i, constraint) in enumerate(block.constraints)
        g[i] = eval_function(constraint, x, block.parameters)
    end
    return
end

function MOI.jacobian_structure(block::QPBlockData)
    J = Tuple{Int,Int}[]
    for (row, constraint) in enumerate(block.constraints)
        append_sparse_gradient_structure!(constraint, J, row)
    end
    return J
end

function MOI.eval_constraint_jacobian(
    block::QPBlockData{T},
    J::AbstractVector{T},
    x::AbstractVector{T},
) where {T}
    i = 1
    for constraint in block.constraints
        ∇f = view(J, i:length(J))
        i += eval_sparse_gradient(∇f, constraint, x, block.parameters)
    end
    return i
end

function MOI.hessian_lagrangian_structure(block::QPBlockData)
    H = Tuple{Int,Int}[]
    append_sparse_hessian_structure!(block.objective, H)
    for constraint in block.constraints
        append_sparse_hessian_structure!(constraint, H)
    end
    return H
end

function MOI.eval_hessian_lagrangian(
    block::QPBlockData{T},
    H::AbstractVector{T},
    x::AbstractVector{T},
    σ::T,
    μ::AbstractVector{T},
) where {T}
    i = 1
    i += eval_sparse_hessian(H, block.objective, σ)
    for (row, constraint) in enumerate(block.constraints)
        ∇²f = view(H, i:length(H))
        i += eval_sparse_hessian(∇²f, constraint, μ[row])
    end
    return i
end