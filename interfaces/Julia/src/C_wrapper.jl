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

function uno_objective(number_variables::Cint, x::Ptr{Float64}, objective_value::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  if isnothing(_user_data.user_model)
    _f = _user_data.eval_objective(_x)::Float64
  else
    _f = _user_data.eval_objective(_user_data.user_model, _x)::Float64
  end
  unsafe_store!(objective_value, _f)
  return Cint(0)
end

function uno_constraints(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, constraint_values::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _c = unsafe_wrap(Array, constraint_values, number_constraints)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  if isnothing(_user_data.user_model)
    _user_data.eval_constraints(_c, _x)
  else
    _user_data.eval_constraints(_user_data.user_model, _c, _x)
  end
  return Cint(0)
end

function uno_objective_gradient(number_variables::Cint, x::Ptr{Float64}, gradient::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _g = unsafe_wrap(Array, gradient, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  if isnothing(_user_data.user_model)
    _user_data.eval_gradient(_g, _x)
  else
    _user_data.eval_gradient(_user_data.user_model, _g, _x)
  end
  return Cint(0)
end

function uno_jacobian(number_variables::Cint, number_jacobian_nonzeros::Cint, x::Ptr{Float64}, jacobian::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _jvals = unsafe_wrap(Array, jacobian, number_jacobian_nonzeros)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  if isnothing(_user_data.user_model)
    _user_data.eval_jacobian(_jvals, _x)
  else
    _user_data.eval_jacobian(_user_data.user_model, _jvals, _x)
  end
  return Cint(0)
end

function uno_lagrangian_hessian(number_variables::Cint, number_constraints::Cint, number_hessian_nonzeros::Cint, x::Ptr{Float64}, objective_multiplier::Float64, multipliers::Ptr{Float64}, hessian::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _multipliers = unsafe_wrap(Array, hessian, number_constraints)
  _hvals = unsafe_wrap(Array, hessian, number_hessian_nonzeros)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  if isnothing(_user_data.user_model)
    _user_data.eval_hessian(_hvals, _x, _multipliers, objective_multiplier)
  else
    _user_data.eval_hessian(_user_data.user_model, _hvals, _x, _multipliers, objective_multiplier)
  end
  return Cint(0)
end

function uno_jacobian_operator(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, evaluate_at_x::Bool, vector::Ptr{Float64}, result::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _v = unsafe_wrap(Array, vector, number_variables)
  _Jv = unsafe_wrap(Array, result, number_constraints)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  if isnothing(_user_data.user_model)
    _user_data.eval_Jv(_Jv, _x, _v, evaluate_at_x)
  else
    _user_data.eval_Jv(_user_data.user_model, _Jv, _x, _v, evaluate_at_x)
  end
  return Cint(0)
end

function uno_jacobian_transposed_operator(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, evaluate_at_x::Bool, vector::Ptr{Float64}, result::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _v = unsafe_wrap(Array, vector, number_constraints)
  _Jtv = unsafe_wrap(Array, result, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  if isnothing(_user_data.user_model)
    _user_data.eval_Jtv(_Jtv, _x, _v, evaluate_at_x)
  else
    _user_data.eval_Jtv(_user_data.user_model, _Jtv, _x, _v, evaluate_at_x)
  end
  return Cint(0)
end

function uno_lagrangian_hessian_operator(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, evaluate_at_x::Bool, objective_multiplier::Float64, multipliers::Ptr{Float64}, vector::Ptr{Float64}, result::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _multipliers = unsafe_wrap(Array, multipliers, number_constraints)
  _v = unsafe_wrap(Array, vector, number_variables)
  _Hv = unsafe_wrap(Array, result, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  if isnothing(_user_data.user_model)
    _user_data.eval_Hv(_Hv, _x, objective_multiplier, _multipliers, _v, evaluate_at_x)
  else
    _user_data.eval_Hv(_user_data.user_model, _Hv, _x, objective_multiplier, _multipliers, _v, evaluate_at_x)
  end
  return Cint(0)
end

mutable struct UnoModel{M}
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
  eval_hessian::Function
  eval_Jv::Union{Function,Nothing}
  eval_Jtv::Union{Function,Nothing}
  eval_Hv::Union{Function,Nothing}
  # User data
  user_model::M
end

function uno_destroy_model(model::UnoModel)
  if model.c_model != C_NULL
    uno_destroy_model(model.c_model)
    model.c_model = C_NULL
  end
end

Base.unsafe_convert(::Type{Ptr{Cvoid}}, model::UnoModel) = model.c_model

function uno_model(
  problem_type::Char,
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
  eval_Hv::Union{Function,Nothing};
  hessian_triangle::Char='L',
  lagrangian_sign::Float64=1.0,
  user_model=nothing
)
  @assert nvar == length(lvar) == length(uvar)
  @assert ncon == length(lcon) == length(ucon)

  # 'L' for linear, 'Q' for quadratic, 'N' for nonlinear
  @assert problem_type == 'L' || problem_type == 'Q' || problem_type == 'N'
  optimization_sense = minimize ? Cint(1) : Cint(-1)

  base_indexing = Cint(1)  # Fortran-style indexing
  c_model = uno_create_model(problem_type, Cint(nvar), lvar, uvar, base_indexing)
  (c_model == C_NULL) && error("Failed to construct Uno model for some unknown reason.")
  model = UnoModel(c_model, nvar, ncon, eval_objective, eval_constraints, eval_gradient,
                   eval_jacobian, eval_hessian, eval_Jv, eval_Jtv, eval_Hv, user_model)

  user_data = pointer_from_objref(model)::Ptr{Cvoid}
  flag = uno_set_user_data(c_model, user_data)
  flag || error("Failed to set user data via uno_set_user_data.")

  eval_objective_c = @cfunction(uno_objective, Cint, (Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
  eval_gradient_c = @cfunction(uno_objective_gradient, Cint, (Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
  flag = uno_set_objective(c_model, optimization_sense, eval_objective_c, eval_gradient_c)
  flag || error("Failed to set objective and gradient via uno_set_objective.")

  if nnzh > 0
    eval_hessian_c = @cfunction(uno_lagrangian_hessian, Cint, (Cint, Cint, Cint, Ptr{Float64}, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
    flag = uno_set_lagrangian_hessian(c_model, Cint(nnzh), hessian_triangle, hrows, hcols, eval_hessian_c, lagrangian_sign)
    flag || error("Failed to set Lagrangian Hessian via uno_set_lagrangian_hessian.")

    if !isnothing(eval_Hv)
      eval_Hv_c = @cfunction(uno_lagrangian_hessian_operator, Cint, (Cint, Cint, Ptr{Float64}, Bool, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
      flag = uno_set_lagrangian_hessian_operator(c_model, Cint(nnzh), eval_Hv_c, lagrangian_sign)
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

  finalizer(uno_destroy_model, model)
  return model
end


function uno_set_initial_primal_iterate(model::UnoModel, initial_primal_iterate::Vector{Float64})
  @assert model.nvar == length(initial_primal_iterate)
  flag = uno_set_initial_primal_iterate(model.c_model, initial_primal_iterate)
  flag || error("Failed to set initial primal iterate via uno_set_initial_primal_iterate.")
  return
end

function uno_set_initial_dual_iterate(model::UnoModel, initial_dual_iterate::Vector{Float64})
  @assert model.ncon == length(initial_dual_iterate)
  flag = uno_set_initial_dual_iterate(model.c_model, initial_dual_iterate)
  flag || error("Failed to set initial dual iterate via uno_set_initial_dual_iterate.")
  return
end

mutable struct UnoSolver
  # Reference to the internal C solver of Uno
  c_solver::Ptr{Cvoid}
end

function uno_destroy_solver(solver::UnoSolver)
  if solver.c_solver != C_NULL
    uno_destroy_solver(solver.c_solver)
    solver.c_solver = C_NULL
  end
end

Base.unsafe_convert(::Type{Ptr{Cvoid}}, solver::UnoSolver) = solver.c_solver

function uno_solver(preset::String; kwargs...)
  c_solver = uno_create_solver()
  (c_solver == C_NULL) && error("Failed to construct Uno solver for some unknown reason.")
  solver = UnoSolver(c_solver)

  # pass options to Uno
  uno_set_solver_preset(c_solver, preset)
  for (k, v) in kwargs
    if typeof(v) <: String
      uno_set_solver_option(c_solver, string(k), v)
    else
      @warn "$k does not seem to be a valid Uno option."
    end
  end

  finalizer(uno_destroy_solver, solver)
  return solver
end

function uno_optimize(solver::UnoSolver, model::UnoModel)
  @assert (solver.c_solver != C_NULL) && (model.c_model != C_NULL)
  uno_optimize(solver.c_solver, model.c_model)
end

function uno_get_cpu_time(solver::UnoSolver)
  @assert solver.c_solver != C_NULL
  uno_get_cpu_time(solver.c_solver)
end

function uno_get_number_iterations(solver::UnoSolver)
  @assert solver.c_solver != C_NULL
  uno_get_number_iterations(solver.c_solver)
end

function uno_get_number_subproblem_solved_evaluations(solver::UnoSolver)
  @assert solver.c_solver != C_NULL
  uno_get_number_subproblem_solved_evaluations(solver.c_solver)
end

function uno_set_solver_option(solver::UnoSolver, option_name::String, option_value::String)
  uno_set_solver_option(solver.c_solver, option_name, option_value)
end

function uno_set_solver_preset(solver::UnoSolver, preset_name::String)
  uno_set_solver_preset(solver.c_solver, preset_name)
end

function uno_get_optimization_status(solver::UnoSolver)
  return uno_get_optimization_status(solver.c_solver)
end

function uno_get_solution_status(solver::UnoSolver)
  return uno_get_solution_status(solver.c_solver)
end

function uno_get_solution_objective(solver::UnoSolver)
  return uno_get_solution_objective(solver.c_solver)
end

function uno_get_primal_solution(solver::UnoSolver, primal_solution::Vector{Float64})
  uno_get_primal_solution(solver.c_solver, primal_solution)
  return primal_solution
end

function uno_get_constraint_dual_solution(solver::UnoSolver, constraint_dual_solution::Vector{Float64})
  uno_get_constraint_dual_solution(solver.c_solver, constraint_dual_solution)
  return constraint_dual_solution
end

function uno_get_lower_bound_dual_solution(solver::UnoSolver, lower_bound_dual_solution::Vector{Float64})
  uno_get_lower_bound_dual_solution(solver.c_solver, lower_bound_dual_solution)
  return lower_bound_dual_solution
end

function uno_get_upper_bound_dual_solution(solver::UnoSolver, upper_bound_dual_solution::Vector{Float64})
  uno_get_upper_bound_dual_solution(solver.c_solver, upper_bound_dual_solution)
  return upper_bound_dual_solution
end

function uno_get_solution_primal_feasibility(solver::UnoSolver)
  return uno_get_solution_primal_feasibility(solver.c_solver)
end

function uno_get_solution_dual_feasibility(solver::UnoSolver)
  return uno_get_solution_dual_feasibility(solver.c_solver)
end

function uno_get_solution_complementarity(solver::UnoSolver)
  return uno_get_solution_complementarity(solver.c_solver)
end
