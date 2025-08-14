function uno_get_version(major, minor, patch)
    @ccall libuno.uno_get_version(major::Ptr{Cint}, minor::Ptr{Cint},
                                  patch::Ptr{Cint})::Cvoid
end

function uno_create_model(problem_type, number_variables, variables_lower_bounds,
                          variables_upper_bounds, vector_indexing)
    @ccall libuno.uno_create_model(problem_type::Cchar, number_variables::Cint,
                                   variables_lower_bounds::Ptr{Cdouble},
                                   variables_upper_bounds::Ptr{Cdouble},
                                   vector_indexing::Cint)::Ptr{Cvoid}
end

function uno_set_objective(model, objective_function, objective_gradient)
    @ccall libuno.uno_set_objective(model::Ptr{Cvoid}, objective_function::Ptr{Cvoid},
                                    objective_gradient::Ptr{Cvoid})::Cvoid
end

function uno_set_constraints(model, number_constraints, constraint_functions,
                             constraints_lower_bounds, constraints_upper_bounds,
                             number_jacobian_nonzeros, jacobian_sparsity,
                             constraint_jacobian)
    @ccall libuno.uno_set_constraints(model::Ptr{Cvoid}, number_constraints::Cint,
                                      constraint_functions::Ptr{Cvoid},
                                      constraints_lower_bounds::Ptr{Cdouble},
                                      constraints_upper_bounds::Ptr{Cdouble},
                                      number_jacobian_nonzeros::Cint,
                                      jacobian_sparsity::Ptr{Cvoid},
                                      constraint_jacobian::Ptr{Cvoid})::Cvoid
end

function uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_sparsity,
                                    lagrangian_hessian)
    @ccall libuno.uno_set_lagrangian_hessian(model::Ptr{Cvoid},
                                             number_hessian_nonzeros::Cint,
                                             hessian_sparsity::Ptr{Cvoid},
                                             lagrangian_hessian::Ptr{Cvoid})::Cvoid
end

function uno_set_user_data(model, user_data)
    @ccall libuno.uno_set_user_data(model::Ptr{Cvoid}, user_data::Ptr{Cvoid})::Cvoid
end

function uno_destroy_model(model)
    @ccall libuno.uno_destroy_model(model::Ptr{Cvoid})::Cvoid
end
