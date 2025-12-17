# Copyright (c) 2025: Alexis Montoison, Charlie Vanaret, other contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE file or at https://opensource.org/licenses/MIT.

function uno_version()
    major = Ref{Cint}()
    minor = Ref{Cint}()
    patch = Ref{Cint}()
    uno_get_version(major, minor, patch)
    return VersionNumber(major[], minor[], patch[])
end

mutable struct Model{M}
  # Reference to the internal C model of Uno
  c_model::Ptr{Cvoid}
  # Number of variables
  nvar::Int
  # Number of constraints
  ncon::Int
  # Callbacks
  eval_objective::Function
  eval_constraints::Function
  eval_gradient::Function
  eval_jacobian::Function
  eval_hessian::Union{Function,Nothing}
  eval_Jv::Union{Function,Nothing}
  eval_Jtv::Union{Function,Nothing}
  eval_Hv::Union{Function,Nothing}
  # User data
  user_model::M
end

function uno_objective(number_variables::Cint, x::Ptr{Float64}, objective_value::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::Model
  _f = _user_data.eval_objective(_user_data.user_model, _x)::Float64
  unsafe_store!(objective_value, _f)
  return Cint(0)
end

function uno_constraints(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, constraint_values::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _c = unsafe_wrap(Array, constraint_values, number_constraints)
  _user_data = unsafe_pointer_to_objref(user_data)::Model
  _user_data.eval_constraints(_user_data.user_model, _c, _x)
  return Cint(0)
end

function uno_objective_gradient(number_variables::Cint, x::Ptr{Float64}, gradient::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _g = unsafe_wrap(Array, gradient, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::Model
  _user_data.eval_gradient(_user_data.user_model, _g, _x)
  return Cint(0)
end

function uno_jacobian(number_variables::Cint, number_jacobian_nonzeros::Cint, x::Ptr{Float64}, jacobian::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _jvals = unsafe_wrap(Array, jacobian, number_jacobian_nonzeros)
  _user_data = unsafe_pointer_to_objref(user_data)::Model
  _user_data.eval_jacobian(_user_data.user_model, _jvals, _x)
  return Cint(0)
end

function uno_lagrangian_hessian(number_variables::Cint, number_constraints::Cint, number_hessian_nonzeros::Cint, x::Ptr{Float64}, objective_multiplier::Float64, multipliers::Ptr{Float64}, hessian::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _multipliers = unsafe_wrap(Array, multipliers, number_constraints)
  _hvals = unsafe_wrap(Array, hessian, number_hessian_nonzeros)
  _user_data = unsafe_pointer_to_objref(user_data)::Model
  _user_data.eval_hessian(_user_data.user_model, _hvals, _x, _multipliers, objective_multiplier)
  return Cint(0)
end

function uno_jacobian_operator(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, evaluate_at_x::Bool, vector::Ptr{Float64}, result::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _v = unsafe_wrap(Array, vector, number_variables)
  _Jv = unsafe_wrap(Array, result, number_constraints)
  _user_data = unsafe_pointer_to_objref(user_data)::Model
  _user_data.eval_Jv(_user_data.user_model, _Jv, _x, _v, evaluate_at_x)
  return Cint(0)
end

function uno_jacobian_transposed_operator(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, evaluate_at_x::Bool, vector::Ptr{Float64}, result::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _v = unsafe_wrap(Array, vector, number_constraints)
  _Jtv = unsafe_wrap(Array, result, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::Model
  _user_data.eval_Jtv(_user_data.user_model, _Jtv, _x, _v, evaluate_at_x)
  return Cint(0)
end

function uno_lagrangian_hessian_operator(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, evaluate_at_x::Bool, objective_multiplier::Float64, multipliers::Ptr{Float64}, vector::Ptr{Float64}, result::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _multipliers = unsafe_wrap(Array, multipliers, number_constraints)
  _v = unsafe_wrap(Array, vector, number_variables)
  _Hv = unsafe_wrap(Array, result, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::Model
  _user_data.eval_Hv(_user_data.user_model, _Hv, _x, objective_multiplier, _multipliers, _v, evaluate_at_x)
  return Cint(0)
end

function uno_destroy_model(model::Model)
  if model.c_model != C_NULL
    uno_destroy_model(model.c_model)
    model.c_model = C_NULL
  end
end

Base.unsafe_convert(::Type{Ptr{Cvoid}}, model::Model) = model.c_model

function uno_model(
  problem_type::String,
  minimize::Bool,
  nvar::Int,
  ncon::Int,
  lvar::Vector{Float64},
  uvar::Vector{Float64},
  lcon::Vector{Float64},
  ucon::Vector{Float64},
  jrows::Vector{Cint},
  jcols::Vector{Cint},
  nnzj::Int,
  hrows::Vector{Cint},
  hcols::Vector{Cint},
  nnzh::Int,
  eval_objective::Function,
  eval_constraints::Function,
  eval_gradient::Function,
  eval_jacobian::Function,
  eval_hessian::Union{Function,Nothing},
  eval_Jv::Union{Function,Nothing},
  eval_Jtv::Union{Function,Nothing},
  eval_Hv::Union{Function,Nothing},
  user_model=nothing,
  hessian_triangle::Char='L',
  lagrangian_sign::Float64=1.0,
  x0::Union{Vector{Float64},Nothing}=nothing,
  y0::Union{Vector{Float64},Nothing}=nothing,
)
  @assert nvar == length(lvar) == length(uvar)
  @assert ncon == length(lcon) == length(ucon)
  @assert isnothing(x0) || nvar == length(x0)
  @assert isnothing(y0) || ncon == length(y0)

  # "LP" for linear problem, "QP" for quadratic problem, "NLP" for nonlinear problem
  @assert problem_type == "LP" || problem_type == "QP" || problem_type == "NLP"
  optimization_sense = minimize ? Cint(1) : Cint(-1)

  base_indexing = Cint(1)  # Fortran-style indexing
  c_model = uno_create_model(problem_type, Cint(nvar), lvar, uvar, base_indexing)
  (c_model == C_NULL) && error("Failed to construct Uno model for some unknown reason.")
  model = Model(c_model, nvar, ncon, eval_objective, eval_constraints, eval_gradient,
                eval_jacobian, eval_hessian, eval_Jv, eval_Jtv, eval_Hv, user_model)

  user_data = pointer_from_objref(model)::Ptr{Cvoid}
  flag = uno_set_user_data(c_model, user_data)
  flag || error("Failed to set user data via uno_set_user_data.")

  eval_objective_c = @cfunction(uno_objective, Cint, (Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
  eval_gradient_c = @cfunction(uno_objective_gradient, Cint, (Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
  flag = uno_set_objective(c_model, optimization_sense, eval_objective_c, eval_gradient_c)
  flag || error("Failed to set objective and gradient via uno_set_objective.")

  if nnzh > 0
    if !isnothing(eval_hessian)
      eval_hessian_c = @cfunction(uno_lagrangian_hessian, Cint, (Cint, Cint, Cint, Ptr{Float64}, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
      flag = uno_set_lagrangian_hessian(c_model, Cint(nnzh), hessian_triangle, hrows, hcols, eval_hessian_c, lagrangian_sign)
      flag || error("Failed to set Lagrangian Hessian via uno_set_lagrangian_hessian.")
    end

    if !isnothing(eval_Hv)
      eval_Hv_c = @cfunction(uno_lagrangian_hessian_operator, Cint, (Cint, Cint, Ptr{Float64}, Bool, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
      flag = uno_set_lagrangian_hessian_operator(c_model, eval_Hv_c, lagrangian_sign)
      flag || error("Failed to set Hessian operator via uno_set_lagrangian_hessian_operator.")
    end
  end

  if ncon > 0
    eval_constraints_c = @cfunction(uno_constraints, Cint, (Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
    eval_jacobian_c = @cfunction(uno_jacobian, Cint, (Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
    flag = uno_set_constraints(c_model, Cint(ncon), eval_constraints_c, lcon, ucon, Cint(nnzj), jrows, jcols, eval_jacobian_c)
    flag || error("Failed to set constraints and Jacobian via uno_set_constraints.")

    if !isnothing(eval_Jv)
      eval_Jv_c = @cfunction(uno_jacobian_operator, Cint, (Cint, Cint, Ptr{Float64}, Bool, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
      flag = uno_set_jacobian_operator(c_model, eval_Jv_c)
      flag || error("Failed to set Jacobian operator via uno_set_jacobian_operator.")
    end

    if !isnothing(eval_Jtv)
      eval_Jtv_c = @cfunction(uno_jacobian_transposed_operator, Cint, (Cint, Cint, Ptr{Float64}, Bool, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
      flag = uno_set_jacobian_transposed_operator(c_model, eval_Jtv_c)
      flag || error("Failed to set transposed Jacobian operator via uno_set_jacobian_transposed_operator.")
    end
  end

  if !isnothing(x0)
    uno_set_initial_primal_iterate(c_model, x0)
  end

  if !isnothing(y0)
    uno_set_initial_dual_iterate(c_model, y0)
  end

  finalizer(uno_destroy_model, model)
  return model
end


function uno(
  problem_type::String,
  minimize::Bool,
  nvar::Int,
  ncon::Int,
  lvar::Vector{Float64},
  uvar::Vector{Float64},
  lcon::Vector{Float64},
  ucon::Vector{Float64},
  jrows::Vector{Cint},
  jcols::Vector{Cint},
  nnzj::Int,
  hrows::Vector{Cint},
  hcols::Vector{Cint},
  nnzh::Int,
  eval_objective::Function,
  eval_constraints::Function,
  eval_gradient::Function,
  eval_jacobian::Function,
  eval_hessian::Function,
  eval_Jv::Union{Function,Nothing},
  eval_Jtv::Union{Function,Nothing},
  eval_Hv::Union{Function,Nothing},
  user_model=nothing,
  hessian_triangle::Char='L',
  lagrangian_sign::Float64=1.0,
  x0::Union{Vector{Float64},Nothing}=nothing,
  y0::Union{Vector{Float64},Nothing}=nothing;
  kwargs...
)
  model = uno_model(problem_type, minimize, nvar, ncon, lvar, uvar, lcon, ucon, jrows,
                    jcols, nnzj, hrows, hcols, nnzh, eval_objective, eval_constraints,
                    eval_gradient, eval_jacobian, eval_hessian, eval_Jv, eval_Jtv,
                    eval_Hv, user_model, hessian_triangle, lagrangian_sign, x0, y0)
  solver = uno_solver(; kwargs...)
  uno_optimize(solver, model)
  return model, solver
end

function uno_set_initial_primal_iterate(model::Model, initial_primal_iterate::Vector{Float64})
  @assert model.nvar == length(initial_primal_iterate)
  GC.@preserve model begin
    flag = uno_set_initial_primal_iterate(model.c_model, initial_primal_iterate)
  end
  flag || error("Failed to set initial primal iterate via uno_set_initial_primal_iterate.")
  return
end

function uno_set_initial_dual_iterate(model::Model, initial_dual_iterate::Vector{Float64})
  @assert model.ncon == length(initial_dual_iterate)
  GC.@preserve model begin
    flag = uno_set_initial_dual_iterate(model.c_model, initial_dual_iterate)
  end
  flag || error("Failed to set initial dual iterate via uno_set_initial_dual_iterate.")
  return
end

mutable struct Solver
  # Reference to the internal C solver of Uno
  c_solver::Ptr{Cvoid}
end

function uno_destroy_solver(solver::Solver)
  if solver.c_solver != C_NULL
    uno_destroy_solver(solver.c_solver)
    solver.c_solver = C_NULL
  end
end

Base.unsafe_convert(::Type{Ptr{Cvoid}}, solver::Solver) = solver.c_solver

function uno_solver(; kwargs...)
  c_solver = uno_create_solver()
  (c_solver == C_NULL) && error("Failed to construct Uno solver for some unknown reason.")
  solver = Solver(c_solver)

  # pass options to Uno
  for (k, v) in kwargs
    if v isa String
      @assert uno_set_solver_string_option(c_solver, string(k), v)
    elseif v isa Float64
      @assert uno_set_solver_double_option(c_solver, string(k), v)
    elseif v isa Bool
      @assert uno_set_solver_bool_option(c_solver, string(k), v)
    elseif v isa Integer
      @assert uno_set_solver_integer_option(c_solver, string(k), Cint(v))
    else
      @warn "$k does not seem to be a valid Uno option."
    end
  end

  finalizer(uno_destroy_solver, solver)
  return solver
end
