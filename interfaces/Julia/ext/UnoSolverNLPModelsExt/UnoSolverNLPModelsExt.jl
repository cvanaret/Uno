module UnoSolverNLPModelsExt

import UnoSolver
import NLPModels
import NLPModels: AbstractNLPModel

function nlpmodels_objective(nlp::AbstractNLPModel{Float64, Vector{Float64}}, x::Vector{Float64})
  f = NLPModels.obj(nlp, x)
  return f
end

function nlpmodels_constraints(nlp::AbstractNLPModel{Float64, Vector{Float64}}, c::Vector{Float64}, x::Vector{Float64})
  NLPModels.cons!(nlp, x, c)
  return c
end

function nlpmodels_objective_gradient(nlp::AbstractNLPModel{Float64, Vector{Float64}}, g::Vector{Float64}, x::Vector{Float64})
  NLPModels.grad!(nlp, x, g)
  return g
end

function nlpmodels_jacobian(nlp::AbstractNLPModel{Float64, Vector{Float64}}, jvals::Vector{Float64}, x::Vector{Float64})
  NLPModels.jac_coord!(nlp, x, jvals)
  return jvals
end

function nlpmodels_lagrangian_hessian(nlp::AbstractNLPModel{Float64, Vector{Float64}}, hvals::Vector{Float64}, x::Vector{Float64},
                                      multipliers::Vector{Float64}, objective_multiplier::Float64)
NLPModels.hess_coord!(nlp, x, multipliers, hvals, obj_weight=objective_multiplier)
  return hvals
end

function nlpmodels_jacobian_operator(nlp::AbstractNLPModel{Float64, Vector{Float64}}, Jv::Vector{Float64}, x::Vector{Float64},
                                     v::Vector{Float64}, evaluate_at_x::Bool)
  NLPModels.jprod!(nlp, x, v, Jv)
  return Jv
end

function nlpmodels_jacobian_transposed_operator(nlp::AbstractNLPModel{Float64, Vector{Float64}}, Jtv::Vector{Float64},
                                                x::Vector{Float64}, v::Vector{Float64}, evaluate_at_x::Bool)
  NLPModels.jtprod!(nlp, x, v, Jtv)
  return Jtv
end

function nlpmodels_lagrangian_hessian_operator(nlp::AbstractNLPModel{Float64, Vector{Float64}}, Hv::Vector{Float64},
                                               x::Vector{Float64}, objective_multiplier::Float64, multipliers::Vector{Float64},
                                               v::Vector{Float64}, evaluate_at_x::Bool)
  NLPModels.hprod!(nlp, x, multipliers, v, Hv; obj_weight=objective_multiplier)
  return Hv
end

function UnoSolver.uno_model(nlp::AbstractNLPModel{Float64, Vector{Float64}}, operators_available::Bool=true)
  jrows, jcols = NLPModels.jac_structure(nlp)
  hrows, hcols = NLPModels.hess_structure(nlp)
  problem_type = nlp.meta.islp ? "LP" : "NLP"
  model = UnoSolver.uno_model(
    problem_type,
    nlp.meta.minimize,
    nlp.meta.nvar,
    nlp.meta.ncon,
    nlp.meta.lvar,
    nlp.meta.uvar,
    nlp.meta.lcon,
    nlp.meta.ucon,
    Cint.(jrows),
    Cint.(jcols),
    nlp.meta.nnzj,
    Cint.(hrows),
    Cint.(hcols),
    nlp.meta.nnzh,
    nlpmodels_objective,
    nlpmodels_constraints,
    nlpmodels_objective_gradient,
    nlpmodels_jacobian,
    nlpmodels_lagrangian_hessian,
    operators_available ? nlpmodels_jacobian_operator : nothing,
    operators_available ? nlpmodels_jacobian_transposed_operator : nothing,
    operators_available ? nlpmodels_lagrangian_hessian_operator : nothing,
    nlp,
    'L',
    1.0,
    nlp.meta.x0,
    nlp.meta.y0
  )
  return model
end

function UnoSolver.uno(nlp::AbstractNLPModel{Float64, Vector{Float64}}, operators_available::Bool=true; kwargs...)
  model = UnoSolver.uno_model(nlp, operators_available)
  solver = UnoSolver.uno_solver(; kwargs...)
  UnoSolver.uno_optimize(solver, model)
  return model, solver
end

end  # module UnoSolverNLPModelsExt
