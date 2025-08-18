# Copyright (c) 2025: Alexis Montoison, Charlie Vanaret, other contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE file or at https://opensource.org/licenses/MIT.

function version()
    major = Ref{Cint}()
    minor = Ref{Cint}()
    patch = Ref{Cint}()
    uno_get_version(major, minor, patch)
    return VersionNumber(major[], minor[], patch[])
end

mutable struct UnoModel
  # Reference to the internal C model of Uno
  c_model::Ptr{Cvoid}
  # Callbacks
  eval_objective::Function
  eval_constraints::Function
  eval_gradient::Function
  eval_jacobian::Function
  eval_hessian::Function
end

Base.unsafe_convert(::Type{Ptr{Cvoid}}, uno_model::UnoModel) = uno_model.c_model

function uno_objective(number_variables::Cint, x::Ptr{Float64}, objective_value::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  _f = _user_data.eval_objective(_x)::Float64
  unsafe_store!(objective_value, _f)
  return Cint(0)
end

function uno_constraints(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, constraint_values::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _c = unsafe_wrap(Array, constraint_values, number_constraints)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  _user_data.eval_constraints(_c, _x)
  return Cint(0)
end

function uno_objective_gradient(number_variables::Cint, x::Ptr{Float64}, gradient::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _g = unsafe_wrap(Array, gradient, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  _user_data.eval_gradient(_g, _x)
  return Cint(0)
end

function uno_jacobian(number_variables::Cint, number_jacobian_nonzeros::Cint, x::Ptr{Float64}, jacobian::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _jvals = unsafe_wrap(Array, jacobian, number_jacobian_nonzeros)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  _user_data.eval_jacobian(_jvals, _x)
  return Cint(0)
end

function uno_hessian(number_variables::Cint, number_constraints::Cint, number_hessian_nonzeros::Cint, x::Ptr{Float64}, objective_multiplier::Float64, multipliers::Ptr{Float64}, hessian::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _multipliers = unsafe_wrap(Array, hessian, number_constraints)
  _hvals = unsafe_wrap(Array, hessian, number_Hessian_nonzeros)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  _user_data.eval_hessian(_hvals, _x, _multipliers, objective_multiplier)
  return Cint(0)
end

function uno_jacobian_vector_product(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, evaluate_at_x::Bool, vector::Ptr{Float64}, result::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _v = unsafe_wrap(Array, vector, number_variables)
  _Jv = unsafe_wrap(Array, result, number_constraints)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  # ...
  return Cint(0)
end

function uno_jacobian_transposed_vector_product(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, evaluate_at_x::Bool, vector::Ptr{Float64}, result::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _v = unsafe_wrap(Array, vector, number_constraints)
  _Jtv = unsafe_wrap(Array, result, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  # ...
  return Cint(0)
end

function uno_hessian_vector_product(number_variables::Cint, number_constraints::Cint, x::Ptr{Float64}, evaluate_at_x::Bool, objective_multiplier::Float64, multipliers::Ptr{Float64}, vector::Ptr{Float64}, result::Ptr{Float64}, user_data::Ptr{Cvoid})
  _x = unsafe_wrap(Array, x, number_variables)
  _multipliers = unsafe_wrap(Array, hessian, number_constraints)
  _v = unsafe_wrap(Array, vector, number_variables)
  _Hv = unsafe_wrap(Array, result, number_variables)
  _user_data = unsafe_pointer_to_objref(user_data)::UnoModel
  # ...
  return Cint(0)
end

function uno(
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
)
  @assert nvar == length(lvar) == length(uvar)
  @assert ncon == length(lcon) == length(ucon)

  problem_type = 'N'  # 'L' for linear, 'Q' for quadratic, 'N' for nonlinear
  vector_indexing = Cint(1)  # Fortran-style indexing
  c_model = uno_create_model(problem_type, Cint(nvar), lvar, uvar, vector_indexing)
  (c_model == C_NULL) && error("Failed to construct Uno model for some unknown reason.")
  uno_model = UnoModel(c_model, eval_objective, eval_constraints, eval_gradient, eval_jacobian, eval_hessian)

  user_data = pointer_from_objref(uno_model)::Ptr{Cvoid}
  uno_set_user_data(c_model, user_data)

  eval_objective_c = @cfunction(uno_objective, Cint, (Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
  eval_gradient_c = @cfunction(uno_objective_gradient, Cint, (Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
  uno_set_objective(c_model, eval_objective_c, eval_gradient_c)

  eval_constraints_c = @cfunction(uno_constraints, Cint, (Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
  eval_jacobian_c = @cfunction(uno_jacobian, Cint, (Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
  uno_set_constraints(c_model, Cint(ncon), eval_constraints_c, lcon, ucon, Cint(nnzj), jrows, jcols, eval_jacobian_c)

  eval_hessian_c = @cfunction(uno_hessian, Cint, (Cint, Cint, Cint, Ptr{Float64}, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
  uno_set_lagrangian_hessian(c_model, Cint(nnzh), 'L', hrows, hcols, eval_hessian_c, 1.0)

  # eval_Jv_c = @cfunction(uno_jacobian_vector_product, Cint, (Cint, Cint, Ptr{Float64}, Bool, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
  # eval_Jtv_c = @cfunction(uno_jacobian_transposed_vector_product, Cint, (Cint, Cint, Ptr{Float64}, Bool, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))
  # eval_Hv_c = @cfunction(uno_hessian_vector_product, Cint, (Cint, Cint, Ptr{Float64}, Bool, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}))

  finalizer(uno_destroy_model, uno_model)
  return uno_model
end
