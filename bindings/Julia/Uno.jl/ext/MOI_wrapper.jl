const _PARAMETER_OFFSET = 0x00f0000000000000

_is_parameter(x::MOI.VariableIndex) = x.value >= _PARAMETER_OFFSET

_is_parameter(term::MOI.ScalarAffineTerm) = _is_parameter(term.variable)

function _is_parameter(term::MOI.ScalarQuadraticTerm)
    return _is_parameter(term.variable_1) || _is_parameter(term.variable_2)
end

"""
    _VectorNonlinearOracle(;
        dimension::Int,
        l::Vector{Float64},
        u::Vector{Float64},
        eval_f::Function,
        jacobian_structure::Vector{Tuple{Int,Int}},
        eval_jacobian::Function,
        hessian_lagrangian_structure::Vector{Tuple{Int,Int}} = Tuple{Int,Int}[],
        eval_hessian_lagrangian::Union{Nothing,Function} = nothing,
    ) <: MOI.AbstractVectorSet

The set:
```math
S = \\{x \\in \\mathbb{R}^{dimension}: l \\le f(x) \\le u\\}
```
where ``f`` is defined by the vectors `l` and `u`, and the callback oracles
`eval_f`, `eval_jacobian`, and `eval_hessian_lagrangian`.

!!! warning
    This set is experimental. We will decide by September 30, 2025, whether to
    convert this into the public `Uno.VectorNonlinearOracle`, move it to
    `MOI.VectorNonlinearOracle`, or remove it completely.

## f

The `eval_f` function must have the signature
```julia
eval_f(ret::AbstractVector, x::AbstractVector)::Nothing
```
which fills ``f(x)`` into the dense vector `ret`.

## Jacobian

The `eval_jacobian` function must have the signature
```julia
eval_jacobian(ret::AbstractVector, x::AbstractVector)::Nothing
```
which fills the sparse Jacobian ``\\nabla f(x)`` into `ret`.

The one-indexed sparsity structure must be provided in the `jacobian_structure`
argument.

## Hessian

The `eval_hessian_lagrangian` function is optional.

If `eval_hessian_lagrangian === nothing`, Uno will use a Hessian approximation
instead of the exact Hessian.

If `eval_hessian_lagrangian` is a function, it must have the signature
```julia
eval_hessian_lagrangian(
    ret::AbstractVector,
    x::AbstractVector,
    μ::AbstractVector,
)::Nothing
```
which fills the sparse Hessian of the Lagrangian ``\\sum \\mu_i \\nabla^2 f_i(x)``
into `ret`.

The one-indexed sparsity structure must be provided in the
`hessian_lagrangian_structure` argument.

## Example

To model the set:
```math
\\begin{align}
0 \\le & x^2           \\le 1
0 \\le & y^2 + z^3 - w \\le 0
\\end{align}
```
do
```jldoctest
julia> import Uno

julia> set = _VectorNonlinearOracle(;
           dimension = 3,
           l = [0.0, 0.0],
           u = [1.0, 0.0],
           eval_f = (ret, x) -> begin
               ret[1] = x[2]^2
               ret[2] = x[3]^2 + x[4]^3 - x[1]
               return
           end,
           jacobian_structure = [(1, 2), (2, 1), (2, 3), (2, 4)],
           eval_jacobian = (ret, x) -> begin
               ret[1] = 2.0 * x[2]
               ret[2] = -1.0
               ret[3] = 2.0 * x[3]
               ret[4] = 3.0 * x[4]^2
               return
           end,
           hessian_lagrangian_structure = [(2, 2), (3, 3), (4, 4)],
           eval_hessian_lagrangian = (ret, x, u) -> begin
               ret[1] = 2.0 * u[1]
               ret[2] = 2.0 * u[2]
               ret[3] = 6.0 * x[4] * u[2]
               return
           end,
       );
```
"""
struct _VectorNonlinearOracle <: MOI.AbstractVectorSet
    input_dimension::Int
    output_dimension::Int
    l::Vector{Float64}
    u::Vector{Float64}
    eval_f::Function
    jacobian_structure::Vector{Tuple{Int,Int}}
    eval_jacobian::Function
    hessian_lagrangian_structure::Vector{Tuple{Int,Int}}
    eval_hessian_lagrangian::Union{Nothing,Function}

    function _VectorNonlinearOracle(;
        dimension::Int,
        l::Vector{Float64},
        u::Vector{Float64},
        eval_f::Function,
        jacobian_structure::Vector{Tuple{Int,Int}},
        eval_jacobian::Function,
        # The hessian_lagrangian is optional.
        hessian_lagrangian_structure::Vector{Tuple{Int,Int}} = Tuple{Int,Int}[],
        eval_hessian_lagrangian::Union{Nothing,Function} = nothing,
    )
        @assert length(l) == length(u)
        return new(
            dimension,
            length(l),
            l,
            u,
            eval_f,
            jacobian_structure,
            eval_jacobian,
            hessian_lagrangian_structure,
            eval_hessian_lagrangian,
        )
    end
end

MOI.dimension(s::_VectorNonlinearOracle) = s.input_dimension

MOI.copy(s::_VectorNonlinearOracle) = s

function Base.show(io::IO, s::_VectorNonlinearOracle)
    println(io, "_VectorNonlinearOracle(;")
    println(io, "    dimension = ", s.input_dimension, ",")
    println(io, "    l = ", s.l, ",")
    println(io, "    u = ", s.u, ",")
    println(io, "    ...,")
    print(io, ")")
    return
end

mutable struct _VectorNonlinearOracleCache
    set::_VectorNonlinearOracle
    x::Vector{Float64}
    eval_f_timer::Float64
    eval_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64

    function _VectorNonlinearOracleCache(set::_VectorNonlinearOracle)
        return new(set, zeros(set.input_dimension), 0.0, 0.0, 0.0)
    end
end

"""
    Optimizer()

Create a new Uno optimizer.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Union{Nothing,Uno.UnoModel}
    name::String
    invalid_model::Bool
    silent::Bool
    options::Dict{String,Any}
    solve_time::Float64
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
    callback::Union{Nothing,Function}
    barrier_iterations::Int
    ad_backend::MOI.Nonlinear.AbstractAutomaticDifferentiation
    vector_nonlinear_oracle_constraints::Vector{
        Tuple{MOI.VectorOfVariables,_VectorNonlinearOracleCache},
    }

    function Optimizer()
        return new(
            nothing,
            "",
            false,
            false,
            Dict{String,Any}(),
            NaN,
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
            nothing,
            0,
            MOI.Nonlinear.SparseReverseMode(),
            Tuple{MOI.VectorOfVariables,_VectorNonlinearOracleCache}[],
        )
    end
end

const _SETS = Union{
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.EqualTo{Float64},
    MOI.Interval{Float64},
}

MOI.get(::Optimizer, ::MOI.SolverVersion) = Uno.version()

### _EmptyNLPEvaluator

struct _EmptyNLPEvaluator <: MOI.AbstractNLPEvaluator end

MOI.features_available(::_EmptyNLPEvaluator) = [:Grad, :Jac, :Hess]
MOI.initialize(::_EmptyNLPEvaluator, ::Any) = nothing
MOI.eval_constraint(::_EmptyNLPEvaluator, g, x) = nothing
MOI.jacobian_structure(::_EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::_EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.eval_constraint_jacobian(::_EmptyNLPEvaluator, J, x) = nothing
MOI.eval_hessian_lagrangian(::_EmptyNLPEvaluator, H, x, σ, μ) = nothing

function MOI.empty!(model::Optimizer)
    model.inner = nothing
    # SKIP: model.name
    model.invalid_model = false
    # SKIP: model.silent
    # SKIP: model.options
    model.solve_time = 0.0
    model.sense = MOI.FEASIBILITY_SENSE
    empty!(model.parameters)
    MOI.empty!(model.variables)
    empty!(model.list_of_variable_indices)
    empty!(model.variable_primal_start)
    empty!(model.mult_x_L)
    empty!(model.mult_x_U)
    model.nlp_data = MOI.NLPBlockData([], _EmptyNLPEvaluator(), false)
    model.nlp_dual_start = nothing
    empty!(model.mult_g_nlp)
    model.qp_data = QPBlockData{Float64}()
    model.nlp_model = nothing
    model.callback = nothing
    model.barrier_iterations = 0
    # SKIP: model.ad_backend
    empty!(model.vector_nonlinear_oracle_constraints)
    return
end

function MOI.is_empty(model::Optimizer)
    return MOI.is_empty(model.variables) &&
           isempty(model.variable_primal_start) &&
           isempty(model.mult_x_L) &&
           isempty(model.mult_x_U) &&
           model.nlp_data.evaluator isa _EmptyNLPEvaluator &&
           model.sense == MOI.FEASIBILITY_SENSE &&
           isempty(model.vector_nonlinear_oracle_constraints)
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(model, src)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Uno"

function MOI.supports_add_constrained_variable(
    ::Optimizer,
    ::Type{MOI.Parameter{Float64}},
)
    return true
end

function _init_nlp_model(model)
    if model.nlp_model === nothing
        if !(model.nlp_data.evaluator isa _EmptyNLPEvaluator)
            error("Cannot mix the new and legacy nonlinear APIs")
        end
        model.nlp_model = MOI.Nonlinear.Model()
    end
    return
end

function MOI.add_constrained_variable(
    model::Optimizer,
    set::MOI.Parameter{Float64},
)
    model.inner = nothing
    _init_nlp_model(model)
    p = MOI.VariableIndex(_PARAMETER_OFFSET + length(model.parameters))
    push!(model.list_of_variable_indices, p)
    model.parameters[p] =
        MOI.Nonlinear.add_parameter(model.nlp_model, set.value)
    ci = MOI.ConstraintIndex{MOI.VariableIndex,typeof(set)}(p.value)
    return p, ci
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Parameter{Float64}},
)
    p = MOI.VariableIndex(ci.value)
    return haskey(model.parameters, p)
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Parameter{Float64}},
    set::MOI.Parameter{Float64},
)
    p = model.parameters[MOI.VariableIndex(ci.value)]
    model.nlp_model[p] = set.value
    return
end

_replace_parameters(model::Optimizer, f) = f

function _replace_parameters(model::Optimizer, f::MOI.VariableIndex)
    if _is_parameter(f)
        return model.parameters[f]
    end
    return f
end

function _replace_parameters(model::Optimizer, f::MOI.ScalarAffineFunction)
    if any(_is_parameter, f.terms)
        g = convert(MOI.ScalarNonlinearFunction, f)
        return _replace_parameters(model, g)
    end
    return f
end

function _replace_parameters(model::Optimizer, f::MOI.ScalarQuadraticFunction)
    if any(_is_parameter, f.affine_terms) ||
       any(_is_parameter, f.quadratic_terms)
        g = convert(MOI.ScalarNonlinearFunction, f)
        return _replace_parameters(model, g)
    end
    return f
end

function _replace_parameters(model::Optimizer, f::MOI.ScalarNonlinearFunction)
    for (i, arg) in enumerate(f.args)
        f.args[i] = _replace_parameters(model, arg)
    end
    return f
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{
        <:Union{
            MOI.VariableIndex,
            MOI.ScalarAffineFunction{Float64},
            MOI.ScalarQuadraticFunction{Float64},
            MOI.ScalarNonlinearFunction,
        },
    },
    ::Type{<:_SETS},
)
    return true
end

### MOI.ListOfConstraintTypesPresent

_add_scalar_nonlinear_constraints(ret, ::Nothing) = nothing

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
    ret = MOI.get(model.variables, attr)
    append!(ret, MOI.get(model.qp_data, attr))
    _add_scalar_nonlinear_constraints(ret, model.nlp_model)
    if !isempty(model.vector_nonlinear_oracle_constraints)
        push!(ret, (MOI.VectorOfVariables, _VectorNonlinearOracle))
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
    MOI.set(model, MOI.RawOptimizerAttribute("max_wall_time"), Float64(value))
    return
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, ::Nothing)
    delete!(model.options, "max_wall_time")
    return
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    return get(model.options, "max_wall_time", nothing)
end

### MOI.RawOptimizerAttribute

MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute) = true

function MOI.set(model::Optimizer, p::MOI.RawOptimizerAttribute, value)
    model.options[p.name] = value
    # No need to reset model.inner because this gets handled in optimize!.
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
    push!(model.mult_x_L, nothing)
    push!(model.mult_x_U, nothing)
    model.inner = nothing
    x = MOI.add_variable(model.variables)
    push!(model.list_of_variable_indices, x)
    return x
end

function MOI.is_valid(model::Optimizer, x::MOI.VariableIndex)
    if _is_parameter(x)
        return haskey(model.parameters, x)
    end
    return MOI.is_valid(model.variables, x)
end

function MOI.get(model::Optimizer, ::MOI.ListOfVariableIndices)
    return model.list_of_variable_indices
end

function MOI.get(model::Optimizer, ::MOI.NumberOfVariables)
    return length(model.list_of_variable_indices)
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    return MOI.is_valid(model.variables, ci)
end

function MOI.get(
    model::Optimizer,
    attr::Union{
        MOI.NumberOfConstraints{MOI.VariableIndex,<:_SETS},
        MOI.ListOfConstraintIndices{MOI.VariableIndex,<:_SETS},
    },
)
    return MOI.get(model.variables, attr)
end

function MOI.get(
    model::Optimizer,
    attr::Union{MOI.ConstraintFunction,MOI.ConstraintSet},
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    return MOI.get(model.variables, attr, c)
end

function MOI.add_constraint(model::Optimizer, x::MOI.VariableIndex, set::_SETS)
    index = MOI.add_constraint(model.variables, x, set)
    model.inner = nothing
    return index
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
    set::S,
) where {S<:_SETS}
    MOI.set(model.variables, MOI.ConstraintSet(), ci, set)
    model.inner = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    MOI.delete(model.variables, ci)
    model.inner = nothing
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
    return MOI.is_valid(model.qp_data, ci)
end

function MOI.add_constraint(
    model::Optimizer,
    func::Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
    set::_SETS,
)
    index = MOI.add_constraint(model.qp_data, func, set)
    model.inner = nothing
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
    return MOI.get(model.qp_data, attr)
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
    return MOI.get(model.qp_data, attr, c)
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
    S<:_SETS,
}
    MOI.set(model.qp_data, MOI.ConstraintSet(), ci, set)
    model.inner = nothing
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
    return MOI.get(model.qp_data, attr, c)
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
    MOI.set(model.qp_data, attr, ci, value)
    # No need to reset model.inner, because this gets handled in optimize!.
    return
end

### ScalarNonlinearFunction

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,<:_SETS},
)
    if model.nlp_model === nothing
        return false
    end
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    return MOI.is_valid(model.nlp_model, index)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ListOfConstraintIndices{F,S},
) where {F<:MOI.ScalarNonlinearFunction,S<:_SETS}
    ret = MOI.ConstraintIndex{F,S}[]
    if model.nlp_model === nothing
        return ret
    end
    for (k, v) in model.nlp_model.constraints
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
    if model.nlp_model === nothing
        return 0
    end
    return count(v.set isa S for v in values(model.nlp_model.constraints))
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.ScalarNonlinearFunction,
    s::_SETS,
)
    _init_nlp_model(model)
    if !isempty(model.parameters)
        _replace_parameters(model, f)
    end
    index = MOI.Nonlinear.add_constraint(model.nlp_model, f, s)
    model.inner = nothing
    return MOI.ConstraintIndex{typeof(f),typeof(s)}(index.value)
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
    _init_nlp_model(model)
    if !isempty(model.parameters)
        _replace_parameters(model, func)
    end
    MOI.Nonlinear.set_objective(model.nlp_model, func)
    model.inner = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,<:_SETS},
)
    MOI.throw_if_not_valid(model, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    return model.nlp_model[index].set
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
    set::S,
) where {S<:_SETS}
    MOI.throw_if_not_valid(model, ci)
    index = MOI.Nonlinear.ConstraintIndex(ci.value)
    func = model.nlp_model[index].expression
    model.nlp_model.constraints[index] = MOI.Nonlinear.Constraint(func, set)
    model.inner = nothing
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
    # No need to reset model.inner, because this gets handled in optimize!.
    return
end

### MOI.VectorOfVariables in _VectorNonlinearOracle

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorOfVariables},
    ::Type{_VectorNonlinearOracle},
)
    return true
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables,_VectorNonlinearOracle},
)
    return 1 <= ci.value <= length(model.vector_nonlinear_oracle_constraints)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ListOfConstraintIndices{F,S},
) where {F<:MOI.VectorOfVariables,S<:_VectorNonlinearOracle}
    n = length(model.vector_nonlinear_oracle_constraints)
    return MOI.ConstraintIndex{F,S}.(1:n)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.NumberOfConstraints{F,S},
) where {F<:MOI.VectorOfVariables,S<:_VectorNonlinearOracle}
    return length(model.vector_nonlinear_oracle_constraints)
end

function MOI.add_constraint(
    model::Optimizer,
    f::F,
    s::S,
) where {F<:MOI.VectorOfVariables,S<:_VectorNonlinearOracle}
    model.inner = nothing
    cache = _VectorNonlinearOracleCache(s)
    push!(model.vector_nonlinear_oracle_constraints, (f, cache))
    n = length(model.vector_nonlinear_oracle_constraints)
    return MOI.ConstraintIndex{F,S}(n)
end

function row(
    model::Optimizer,
    ci::MOI.ConstraintIndex{F,S},
) where {F<:MOI.VectorOfVariables,S<:_VectorNonlinearOracle}
    offset = length(model.qp_data)
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
) where {F<:MOI.VectorOfVariables,S<:_VectorNonlinearOracle}
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    f, _ = model.vector_nonlinear_oracle_constraints[ci.value]
    return MOI.get.(model, MOI.VariablePrimal(attr.result_index), f.variables)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{F,S},
) where {F<:MOI.VectorOfVariables,S<:_VectorNonlinearOracle}
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    sign = -_dual_multiplier(model)
    f, s = model.vector_nonlinear_oracle_constraints[ci.value]
    λ = model.inner.mult_g[row(model, ci)]
    J = Tuple{Int,Int}[]
    _jacobian_structure(J, 0, f, s)
    J_val = zeros(length(J))
    _eval_constraint_jacobian(J_val, 0, model.inner.x, f, s)
    dual = zeros(MOI.dimension(s.set))
    # dual = λ' * J(x)
    col_to_index = Dict(x.value => j for (j, x) in enumerate(f.variables))
    for ((row, col), J_rc) in zip(J, J_val)
        dual[col_to_index[col]] += sign * J_rc * λ[row]
    end
    return dual
end

### UserDefinedFunction

MOI.supports(model::Optimizer, ::MOI.UserDefinedFunction) = true

function MOI.set(model::Optimizer, attr::MOI.UserDefinedFunction, args)
    _init_nlp_model(model)
    MOI.Nonlinear.register_operator(
        model.nlp_model,
        attr.name,
        attr.arity,
        args...,
    )
    return
end

### ListOfSupportedNonlinearOperators

function MOI.get(model::Optimizer, attr::MOI.ListOfSupportedNonlinearOperators)
    _init_nlp_model(model)
    return MOI.get(model.nlp_model, attr)
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
    if _is_parameter(vi)
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
    if _is_parameter(vi)
        throw(MOI.SetAttributeNotAllowed(attr, "Variable is a Parameter"))
    end
    MOI.throw_if_not_valid(model, vi)
    model.variable_primal_start[column(vi)] = value
    # No need to reset model.inner, because this gets handled in optimize!.
    return
end

### MOI.ConstraintDualStart

_dual_start(::Optimizer, ::Nothing, ::Int = 1) = 0.0

function _dual_start(model::Optimizer, value::Real, scale::Int = 1)
    return _dual_multiplier(model) * value * scale
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ConstraintDualStart,
    ::Type{MOI.ConstraintIndex{MOI.VariableIndex,S}},
) where {S<:_SETS}
    return true
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
    value::Union{Real,Nothing},
)
    MOI.throw_if_not_valid(model, ci)
    model.mult_x_L[ci.value] = value
    # No need to reset model.inner, because this gets handled in optimize!.
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
)
    MOI.throw_if_not_valid(model, ci)
    return model.mult_x_L[ci.value]
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
    value::Union{Real,Nothing},
)
    MOI.throw_if_not_valid(model, ci)
    model.mult_x_U[ci.value] = value
    # No need to reset model.inner, because this gets handled in optimize!.
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
)
    MOI.throw_if_not_valid(model, ci)
    return model.mult_x_U[ci.value]
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
    value::Union{Real,Nothing},
) where {S<:Union{MOI.EqualTo{Float64},MOI.Interval{Float64}}}
    MOI.throw_if_not_valid(model, ci)
    if value === nothing
        model.mult_x_L[ci.value] = nothing
        model.mult_x_U[ci.value] = nothing
    elseif value >= 0.0
        model.mult_x_L[ci.value] = value
        model.mult_x_U[ci.value] = 0.0
    else
        model.mult_x_L[ci.value] = 0.0
        model.mult_x_U[ci.value] = value
    end
    # No need to reset model.inner, because this gets handled in optimize!.
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:Union{MOI.EqualTo{Float64},MOI.Interval{Float64}}}
    MOI.throw_if_not_valid(model, ci)
    l = model.mult_x_L[ci.value]
    u = model.mult_x_U[ci.value]
    return (l === u === nothing) ? nothing : (l + u)
end

### MOI.NLPBlockDualStart

MOI.supports(::Optimizer, ::MOI.NLPBlockDualStart) = true

function MOI.set(
    model::Optimizer,
    ::MOI.NLPBlockDualStart,
    values::Union{Nothing,Vector},
)
    model.nlp_dual_start = values
    # No need to reset model.inner, because this gets handled in optimize!.
    return
end

MOI.get(model::Optimizer, ::MOI.NLPBlockDualStart) = model.nlp_dual_start

### MOI.NLPBlock

MOI.supports(::Optimizer, ::MOI.NLPBlock) = true

# This may also be set by `optimize!` and contain the block created from
# ScalarNonlinearFunction
MOI.get(model::Optimizer, ::MOI.NLPBlock) = model.nlp_data

function MOI.set(model::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    if model.nlp_model !== nothing
        error("Cannot mix the new and legacy nonlinear APIs")
    end
    model.nlp_data = nlp_data
    model.inner = nothing
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
    model.inner = nothing
    return
end

MOI.get(model::Optimizer, ::MOI.ObjectiveSense) = model.sense

### ObjectiveFunction

function MOI.get(model::Optimizer, attr::MOI.ObjectiveFunctionType)
    if model.nlp_model !== nothing && model.nlp_model.objective !== nothing
        return MOI.ScalarNonlinearFunction
    end
    return MOI.get(model.qp_data, attr)
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
    return MOI.get(model.qp_data, attr)
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
    MOI.set(model.qp_data, attr, func)
    if model.nlp_model !== nothing
        MOI.Nonlinear.set_objective(model.nlp_model, nothing)
    end
    model.inner = nothing
    return
end

function MOI.eval_objective(model::Optimizer, x)
    # TODO(odow): FEASIBILITY_SENSE could produce confusing solver output if
    # a nonzero objective is set.
    if model.sense == MOI.FEASIBILITY_SENSE
        return 0.0
    elseif model.nlp_data.has_objective
        return MOI.eval_objective(model.nlp_data.evaluator, x)::Float64
    end
    return MOI.eval_objective(model.qp_data, x)
end

function MOI.eval_objective_gradient(model::Optimizer, grad, x)
    if model.sense == MOI.FEASIBILITY_SENSE
        grad .= zero(eltype(grad))
    elseif model.nlp_data.has_objective
        MOI.eval_objective_gradient(model.nlp_data.evaluator, grad, x)
    else
        MOI.eval_objective_gradient(model.qp_data, grad, x)
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

function MOI.eval_constraint(model::Optimizer, g, x)
    MOI.eval_constraint(model.qp_data, g, x)
    offset = length(model.qp_data)
    for (f, s) in model.vector_nonlinear_oracle_constraints
        offset = _eval_constraint(g, offset, x, f, s)
    end
    g_nlp = view(g, (offset+1):length(g))
    MOI.eval_constraint(model.nlp_data.evaluator, g_nlp, x)
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

function MOI.jacobian_structure(model::Optimizer)
    J = MOI.jacobian_structure(model.qp_data)
    offset = length(model.qp_data)
    for (f, s) in model.vector_nonlinear_oracle_constraints
        offset = _jacobian_structure(J, offset, f, s)
    end
    if length(model.nlp_data.constraint_bounds) > 0
        J_nlp = MOI.jacobian_structure(
            model.nlp_data.evaluator,
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

function MOI.eval_constraint_jacobian(model::Optimizer, values, x)
    offset = MOI.eval_constraint_jacobian(model.qp_data, values, x)
    offset -= 1  # .qp_data returns one-indexed offset
    for (f, s) in model.vector_nonlinear_oracle_constraints
        offset = _eval_constraint_jacobian(values, offset, x, f, s)
    end
    nlp_values = view(values, (offset+1):length(values))
    MOI.eval_constraint_jacobian(model.nlp_data.evaluator, nlp_values, x)
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

function MOI.hessian_lagrangian_structure(model::Optimizer)
    H = MOI.hessian_lagrangian_structure(model.qp_data)
    for (f, s) in model.vector_nonlinear_oracle_constraints
        _hessian_lagrangian_structure(H, f, s)
    end
    append!(H, MOI.hessian_lagrangian_structure(model.nlp_data.evaluator))
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

function MOI.eval_hessian_lagrangian(model::Optimizer, H, x, σ, μ)
    offset = MOI.eval_hessian_lagrangian(model.qp_data, H, x, σ, μ)
    offset -= 1  # .qp_data returns one-indexed offset
    μ_offset = length(model.qp_data)
    for (f, s) in model.vector_nonlinear_oracle_constraints
        offset, μ_offset =
            _eval_hessian_lagrangian(H, offset, x, μ, μ_offset, f, s)
    end
    H_nlp = view(H, (offset+1):length(H))
    μ_nlp = view(μ, (μ_offset+1):length(μ))
    MOI.eval_hessian_lagrangian(model.nlp_data.evaluator, H_nlp, x, σ, μ_nlp)
    return
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
    model.ad_backend = backend
    return
end

### MOI.optimize!

function _setup_model(model::Optimizer)
    vars = MOI.get(model.variables, MOI.ListOfVariableIndices())
    if isempty(vars)
        # Don't attempt to create a problem because Uno will error.
        model.invalid_model = true
        return
    end
    if model.nlp_model !== nothing
        model.nlp_data = MOI.NLPBlockData(
            MOI.Nonlinear.Evaluator(model.nlp_model, model.ad_backend, vars),
        )
    end
    has_quadratic_constraints =
        any(isequal(_kFunctionTypeScalarQuadratic), model.qp_data.function_type)
    has_nlp_constraints =
        !isempty(model.nlp_data.constraint_bounds) ||
        !isempty(model.vector_nonlinear_oracle_constraints)
    has_hessian = :Hess in MOI.features_available(model.nlp_data.evaluator)
    for (_, s) in model.vector_nonlinear_oracle_constraints
        if s.set.eval_hessian_lagrangian === nothing
            has_hessian = false
            break
        end
    end
    init_feat = [:Grad]
    if has_hessian
        push!(init_feat, :Hess)
    end
    if has_nlp_constraints
        push!(init_feat, :Jac)
    end
    MOI.initialize(model.nlp_data.evaluator, init_feat)

    jacobian_sparsity = MOI.jacobian_structure(model)
    nnzj = length(jacobian_sparsity)
    jrows = Vector{Cint}(undef, nnzj)
    jcols = Vector{Cint}(undef, nnzj)
    for i in 1:nnzj
        jrows[i], jcols[i] = jacobian_sparsity[i]
    end

    hessian_sparsity = has_hessian ? MOI.hessian_lagrangian_structure(model) : Tuple{Int,Int}[]
    nnzh = length(hessian_sparsity)
    hrows = Vector{Cint}(undef, nnzh)
    hcols = Vector{Cint}(undef, nnzh)
    for i in 1:nnzh
        hrows[i], hcols[i] = hessian_sparsity[i]
    end

    moi_objective(model, x) = MOI.eval_objective(model, x)
    moi_objective_gradient(model, g, x) = MOI.eval_objective_gradient(model, g, x)
    moi_constraints(model, c, x) = MOI.eval_constraint(model, c, x)
    moi_jacobian(model, jvals, x) = MOI.eval_constraint_jacobian(model, jvals, x)
    moi_lagrangian_hessian(model, hvals, x, multipliers, objective_multiplier) = MOI.eval_hessian_lagrangian(model, hvals, x, objective_multiplier, multipliers)

    # To be fixed by Alexis one day...
    moi_jacobian_operator(model, Jv, x, v, evaluate_at_x) = nothing
    moi_jacobian_transposed_operator(model, Jtv, x, v, evaluate_at_x) = nothing
    moi_lagrangian_hessian_operator(model, Hv, x, objective_multiplier, multipliers, v, evaluate_at_x) = nothing

    g_L, g_U = copy(model.qp_data.g_L), copy(model.qp_data.g_U)
    for (_, s) in model.vector_nonlinear_oracle_constraints
        append!(g_L, s.set.l)
        append!(g_U, s.set.u)
    end
    for bound in model.nlp_data.constraint_bounds
        push!(g_L, bound.lower)
        push!(g_U, bound.upper)
    end
    nvar = length(vars)
    ncon = length(g_L)
    x0 = zeros(Float64, nvar)
    y0 = zeros(Float64, ncon)
    model.inner = Uno.uno(
        'N',
        model.sense == MOI.MIN_SENSE,
        nvar,
        ncon,
        model.variables.lower,
        model.variables.upper,
        g_L,
        g_U,
        jrows,
        jcols,
        nnzj,
        hrows,
        hcols,
        nnzh,
        x0,
        y0,
        moi_objective,
        moi_constraints,
        moi_objective_gradient,
        moi_jacobian,
        moi_lagrangian_hessian,
        moi_jacobian_operator,
        moi_jacobian_transposed_operator,
        moi_lagrangian_hessian_operator;
        hessian_triangle='L',
        lagrangian_sign=1.0,
        user_model=model,
    )
    return
end

function copy_parameters(model::Optimizer)
    if model.nlp_model === nothing
        return
    end
    empty!(model.qp_data.parameters)
    for (p, index) in model.parameters
        model.qp_data.parameters[p.value] = model.nlp_model[index]
    end
    return
end

function MOI.optimize!(model::Optimizer)
    start_time = time()
    if model.inner === nothing
        _setup_model(model)
    end
    if model.invalid_model
        return
    end
    copy_parameters(model)
    inner = model.inner::Uno.UnoModel

    # Alexis: I need the help of Charlie for that!
    # The default print level is `5`
    # Uno.AddUnoIntOption(inner, "print_level", model.silent ? 0 : 5)
    Uno.uno_set_solver_option(model.inner, "print_solution", "yes")

    # Other misc options that over-ride the ones set above.
    # for (name, value) in model.options
    #     if value isa String
    #         Uno.AddUnoStrOption(inner, name, value)
    #     elseif value isa Integer
    #         Uno.AddUnoIntOption(inner, name, value)
    #     elseif value isa Float64
    #         Uno.AddUnoNumOption(inner, name, value)
    #     else
    #         error(
    #             "Unable to add option `\"$name\"` with the value " *
    #             "`$value::$(typeof(value))`. The value must be a `::String`, " *
    #             "`::Integer`, or `::Float64`.",
    #         )
    #     end
    # end
    # Initialize the starting point, projecting variables from 0 onto their
    # bounds if VariablePrimalStart is not provided.
    # for i in 1:length(model.variable_primal_start)
    #     inner.x[i] = if model.variable_primal_start[i] !== nothing
    #         model.variable_primal_start[i]
    #     else
    #         clamp(0.0, model.variables.lower[i], model.variables.upper[i])
    #     end
    # end
    # for (i, start) in enumerate(model.qp_data.mult_g)
    #     inner.mult_g[i] = _dual_start(model, start, -1)
    # end
    # offset = length(model.qp_data.mult_g)
    # if model.nlp_dual_start === nothing
    #     inner.mult_g[(offset+1):end] .= 0.0
    #     for (key, val) in model.mult_g_nlp
    #         inner.mult_g[offset+key.value] = val
    #     end
    # else
    #     for (i, start) in enumerate(model.nlp_dual_start::Vector{Float64})
    #         inner.mult_g[offset+i] = _dual_start(model, start, -1)
    #     end
    # end
    # for i in 1:inner.n
    #     inner.mult_x_L[i] = _dual_start(model, model.mult_x_L[i])
    #     inner.mult_x_U[i] = _dual_start(model, model.mult_x_U[i], -1)
    # end
    model.barrier_iterations = 0
    # Clear timers
    for (_, s) in model.vector_nonlinear_oracle_constraints
        s.eval_f_timer = 0.0
        s.eval_jacobian_timer = 0.0
        s.eval_hessian_lagrangian_timer = 0.0
    end
    Uno.uno_optimize(model.inner)
    model.solve_time = time() - start_time
    return
end

const _STATUS_CODES = Dict{
    Uno.ApplicationReturnStatus, Tuple{MOI.TerminationStatusCode, MOI.ResultStatusCode}
}(
    Uno.UNO_SUCCESS           => (MOI.LOCALLY_SOLVED,  MOI.FEASIBLE_POINT       ),
    Uno.UNO_ITERATION_LIMIT   => (MOI.ITERATION_LIMIT, MOI.UNKNOWN_RESULT_STATUS),
    Uno.UNO_TIME_LIMIT        => (MOI.TIME_LIMIT,      MOI.UNKNOWN_RESULT_STATUS),
    Uno.UNO_EVALUATION_ERROR  => (MOI.NUMERICAL_ERROR, MOI.UNKNOWN_RESULT_STATUS),
    Uno.UNO_ALGORITHMIC_ERROR => (MOI.OTHER_ERROR,     MOI.UNKNOWN_RESULT_STATUS),
)

### MOI.ResultCount

# Uno always has an iterate available.
function MOI.get(model::Optimizer, ::MOI.ResultCount)
    return (model.inner !== nothing) ? 1 : 0
end

### MOI.TerminationStatus

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.invalid_model
        return MOI.INVALID_MODEL
    elseif model.inner === nothing
        return MOI.OPTIMIZE_NOT_CALLED
    end
    status, _ = _STATUS_CODES[Uno.ApplicationReturnStatus(model.inner.status)]
    return status
end

### MOI.RawStatusString

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    if model.invalid_model
        return "The model has no variable"
    elseif model.inner === nothing
        return "Optimize not called"
    end
    return string(Uno.ApplicationReturnStatus(model.inner.status))
end

### MOI.PrimalStatus

function _manually_evaluated_primal_status(model::Optimizer)
    x, g = model.inner.x, model.inner.g
    m, n = length(g), length(x)
    x_L, x_U = model.variables.lower, model.variables.upper
    g_L, g_U = copy(model.qp_data.g_L), copy(model.qp_data.g_U)
    # Assuming constraints are guaranteed to be in the order [qp_cons, nlp_cons]
    for bound in model.nlp_data.constraint_bounds
        push!(g_L, bound.lower)
        push!(g_U, bound.upper)
    end
    # 1e-8 is the default tolerance
    tol = get(model.options, "tol", 1e-8)
    if all(x_L[i] - tol <= x[i] <= x_U[i] + tol for i in 1:n) &&
       all(g_L[i] - tol <= g[i] <= g_U[i] + tol for i in 1:m)
        return MOI.FEASIBLE_POINT
    end
    # 1e-6 is the default acceptable tolerance
    atol = get(model.options, "acceptable_tol", 1e-6)
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
    _, status = _STATUS_CODES[Uno.ApplicationReturnStatus(model.inner.status)]
    if status == MOI.UNKNOWN_RESULT_STATUS
        # Not sure. RestorationFailure can terminate at a feasible (but
        # non-stationary) point.
        return _manually_evaluated_primal_status(model)
    end
    return status
end

### MOI.DualStatus

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    if !(1 <= attr.result_index <= MOI.get(model, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end
    _, status = _STATUS_CODES[Uno.ApplicationReturnStatus(model.inner.status)]
    return status
end

### MOI.SolveTimeSec

MOI.get(model::Optimizer, ::MOI.SolveTimeSec) = model.solve_time

### MOI.BarrierIterations

MOI.get(model::Optimizer, ::MOI.BarrierIterations) = model.barrier_iterations

### MOI.ObjectiveValue

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return model.inner.obj_val
end

### MOI.VariablePrimal

function MOI.get(
    model::Optimizer,
    attr::MOI.VariablePrimal,
    vi::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, vi)
    if _is_parameter(vi)
        p = model.parameters[vi]
        return model.nlp_model[p]
    end
    return model.inner.x[column(vi)]
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
    offset = length(model.qp_data)
    for (_, s) in model.vector_nonlinear_oracle_constraints
        offset += s.set.output_dimension
    end
    return offset + ci.value
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{F,<:_SETS},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
        MOI.ScalarNonlinearFunction,
    },
}
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    return model.inner.g[row(model, ci)]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:_SETS},
)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    return model.inner.x[ci.value]
end

### MOI.ConstraintDual

_dual_multiplier(model::Optimizer) = model.sense == MOI.MIN_SENSE ? 1.0 : -1.0

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{F,<:_SETS},
) where {
    F<:Union{
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
        MOI.ScalarNonlinearFunction,
    },
}
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    s = -_dual_multiplier(model)
    return s * model.inner.mult_g[row(model, ci)]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    rc = model.inner.mult_x_L[ci.value] - model.inner.mult_x_U[ci.value]
    return min(0.0, _dual_multiplier(model) * rc)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    rc = model.inner.mult_x_L[ci.value] - model.inner.mult_x_U[ci.value]
    return max(0.0, _dual_multiplier(model) * rc)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}},
)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    rc = model.inner.mult_x_L[ci.value] - model.inner.mult_x_U[ci.value]
    return _dual_multiplier(model) * rc
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Interval{Float64}},
)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    rc = model.inner.mult_x_L[ci.value] - model.inner.mult_x_U[ci.value]
    return _dual_multiplier(model) * rc
end

### MOI.NLPBlockDual

function MOI.get(model::Optimizer, attr::MOI.NLPBlockDual)
    MOI.check_result_index_bounds(model, attr)
    s = -_dual_multiplier(model)
    return s .* model.inner.mult_g[(length(model.qp_data)+1):end]
end

### We need to Oscar why we need this routine

function MOI.get(
    model::Optimizer,
    ::MOI.CallbackVariablePrimal,
    x::MOI.VariableIndex,
)
    return model.inner.x[column(x)]
end