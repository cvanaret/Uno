function uno_get_version(major, minor, patch)
    @ccall libuno.uno_get_version(major::Ptr{Int32}, minor::Ptr{Int32},
                                  patch::Ptr{Int32})::Cvoid
end

function uno_create_model(problem_type, number_variables, variables_lower_bounds,
                          variables_upper_bounds, vector_indexing)
    @ccall libuno.uno_create_model(problem_type::Cchar, number_variables::Int32,
                                   variables_lower_bounds::Ptr{Cdouble},
                                   variables_upper_bounds::Ptr{Cdouble},
                                   vector_indexing::Int32)::Ptr{Cvoid}
end

function uno_set_objective(model, objective_sense, objective_function, objective_gradient)
    @ccall libuno.uno_set_objective(model::Ptr{Cvoid}, objective_sense::Cdouble,
                                    objective_function::Ptr{Cvoid},
                                    objective_gradient::Ptr{Cvoid})::Cvoid
end

function uno_set_constraints(model, number_constraints, constraint_functions,
                             constraints_lower_bounds, constraints_upper_bounds,
                             number_jacobian_nonzeros, jacobian_row_indices,
                             jacobian_column_indices, constraint_jacobian)
    @ccall libuno.uno_set_constraints(model::Ptr{Cvoid}, number_constraints::Int32,
                                      constraint_functions::Ptr{Cvoid},
                                      constraints_lower_bounds::Ptr{Cdouble},
                                      constraints_upper_bounds::Ptr{Cdouble},
                                      number_jacobian_nonzeros::Int32,
                                      jacobian_row_indices::Ptr{Int32},
                                      jacobian_column_indices::Ptr{Int32},
                                      constraint_jacobian::Ptr{Cvoid})::Cvoid
end

function uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part,
                                    hessian_row_indices, hessian_column_indices,
                                    lagrangian_hessian, lagrangian_sign_convention)
    @ccall libuno.uno_set_lagrangian_hessian(model::Ptr{Cvoid},
                                             number_hessian_nonzeros::Int32,
                                             hessian_triangular_part::Cchar,
                                             hessian_row_indices::Ptr{Int32},
                                             hessian_column_indices::Ptr{Int32},
                                             lagrangian_hessian::Ptr{Cvoid},
                                             lagrangian_sign_convention::Cdouble)::Cvoid
end

function uno_set_user_data(model, user_data)
    @ccall libuno.uno_set_user_data(model::Ptr{Cvoid}, user_data::Ptr{Cvoid})::Cvoid
end

function uno_destroy_model(model)
    @ccall libuno.uno_destroy_model(model::Ptr{Cvoid})::Cvoid
end
