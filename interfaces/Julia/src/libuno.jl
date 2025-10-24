const UNO_MINIMIZE = Cint(1)
const UNO_MAXIMIZE = Cint(-1)

const UNO_MULTIPLIER_POSITIVE = Cdouble(1.0)
const UNO_MULTIPLIER_NEGATIVE = Cdouble(-1.0)

const UNO_PROBLEM_LINEAR = "LP"
const UNO_PROBLEM_QUADRATIC = "QP"
const UNO_PROBLEM_NONLINEAR = "NLP"

const UNO_ZERO_BASED_INDEXING = Cint(0)
const UNO_ONE_BASED_INDEXING = Cint(1)

const UNO_LOWER_TRIANGLE = Cchar('L')
const UNO_UPPER_TRIANGLE = Cchar('U')

const UNO_OPTION_TYPE_INTEGER = Cint(0)
const UNO_OPTION_TYPE_DOUBLE = Cint(1)
const UNO_OPTION_TYPE_BOOL = Cint(2)
const UNO_OPTION_TYPE_STRING = Cint(3)
const UNO_OPTION_TYPE_NOT_FOUND = Cint(-1)

const UNO_SUCCESS = Cint(0)
const UNO_ITERATION_LIMIT = Cint(1)
const UNO_TIME_LIMIT = Cint(2)
const UNO_EVALUATION_ERROR = Cint(3)
const UNO_ALGORITHMIC_ERROR = Cint(4)
const UNO_USER_TERMINATION = Cint(5)

const UNO_NOT_OPTIMAL = Cint(0)
const UNO_FEASIBLE_KKT_POINT = Cint(1)
const UNO_FEASIBLE_FJ_POINT = Cint(2)
const UNO_INFEASIBLE_STATIONARY_POINT = Cint(3)
const UNO_FEASIBLE_SMALL_STEP = Cint(4)
const UNO_INFEASIBLE_SMALL_STEP = Cint(5)
const UNO_UNBOUNDED = Cint(6)

function uno_get_version(major, minor, patch)
    @ccall libuno.uno_get_version(major::Ptr{Int32}, minor::Ptr{Int32},
                                  patch::Ptr{Int32})::Cvoid
end

function uno_create_model(problem_type, number_variables, variables_lower_bounds,
                          variables_upper_bounds, base_indexing)
    @ccall libuno.uno_create_model(problem_type::Cstring, number_variables::Int32,
                                   variables_lower_bounds::Ptr{Cdouble},
                                   variables_upper_bounds::Ptr{Cdouble},
                                   base_indexing::Int32)::Ptr{Cvoid}
end

function uno_set_objective(model, optimization_sense, objective_function,
                           objective_gradient)
    @ccall libuno.uno_set_objective(model::Ptr{Cvoid}, optimization_sense::Int32,
                                    objective_function::Ptr{Cvoid},
                                    objective_gradient::Ptr{Cvoid})::Bool
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
                                      constraint_jacobian::Ptr{Cvoid})::Bool
end

function uno_set_jacobian_operator(model, jacobian_operator)
    @ccall libuno.uno_set_jacobian_operator(model::Ptr{Cvoid},
                                            jacobian_operator::Ptr{Cvoid})::Bool
end

function uno_set_jacobian_transposed_operator(model, jacobian_transposed_operator)
    @ccall libuno.uno_set_jacobian_transposed_operator(model::Ptr{Cvoid},
                                                       jacobian_transposed_operator::Ptr{Cvoid})::Bool
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
                                             lagrangian_sign_convention::Cdouble)::Bool
end

function uno_set_lagrangian_hessian_operator(model, lagrangian_hessian_operator,
                                             lagrangian_sign_convention)
    @ccall libuno.uno_set_lagrangian_hessian_operator(model::Ptr{Cvoid},
                                                      lagrangian_hessian_operator::Ptr{Cvoid},
                                                      lagrangian_sign_convention::Cdouble)::Bool
end

function uno_set_user_data(model, user_data)
    @ccall libuno.uno_set_user_data(model::Ptr{Cvoid}, user_data::Ptr{Cvoid})::Bool
end

function uno_set_initial_primal_iterate_component(model, index, initial_primal_component)
    @ccall libuno.uno_set_initial_primal_iterate_component(model::Ptr{Cvoid}, index::Int32,
                                                           initial_primal_component::Cdouble)::Bool
end

function uno_set_initial_dual_iterate_component(model, index, initial_dual_component)
    @ccall libuno.uno_set_initial_dual_iterate_component(model::Ptr{Cvoid}, index::Int32,
                                                         initial_dual_component::Cdouble)::Bool
end

function uno_set_initial_primal_iterate(model, initial_primal_iterate)
    @ccall libuno.uno_set_initial_primal_iterate(model::Ptr{Cvoid},
                                                 initial_primal_iterate::Ptr{Cdouble})::Bool
end

function uno_set_initial_dual_iterate(model, initial_dual_iterate)
    @ccall libuno.uno_set_initial_dual_iterate(model::Ptr{Cvoid},
                                               initial_dual_iterate::Ptr{Cdouble})::Bool
end

function uno_create_solver()
    @ccall libuno.uno_create_solver()::Ptr{Cvoid}
end

function uno_set_solver_integer_option(solver, option_name, option_value)
    @ccall libuno.uno_set_solver_integer_option(solver::Ptr{Cvoid}, option_name::Cstring,
                                                option_value::Int32)::Bool
end

function uno_set_solver_double_option(solver, option_name, option_value)
    @ccall libuno.uno_set_solver_double_option(solver::Ptr{Cvoid}, option_name::Cstring,
                                               option_value::Cdouble)::Bool
end

function uno_set_solver_bool_option(solver, option_name, option_value)
    @ccall libuno.uno_set_solver_bool_option(solver::Ptr{Cvoid}, option_name::Cstring,
                                             option_value::Bool)::Bool
end

function uno_set_solver_string_option(solver, option_name, option_value)
    @ccall libuno.uno_set_solver_string_option(solver::Ptr{Cvoid}, option_name::Cstring,
                                               option_value::Cstring)::Bool
end

function uno_get_solver_option_type(solver, option_name)
    @ccall libuno.uno_get_solver_option_type(solver::Ptr{Cvoid},
                                             option_name::Cstring)::Int32
end

function uno_load_solver_option_file(solver, file_name)
    @ccall libuno.uno_load_solver_option_file(solver::Ptr{Cvoid}, file_name::Cstring)::Bool
end

function uno_set_solver_preset(solver, preset_name)
    @ccall libuno.uno_set_solver_preset(solver::Ptr{Cvoid}, preset_name::Cstring)::Bool
end

function uno_set_solver_callbacks(solver, notify_acceptable_iterate_callback,
                                  user_termination_callback, user_data)
    @ccall libuno.uno_set_solver_callbacks(solver::Ptr{Cvoid},
                                           notify_acceptable_iterate_callback::Ptr{Cvoid},
                                           user_termination_callback::Ptr{Cvoid},
                                           user_data::Ptr{Cvoid})::Bool
end

function uno_set_logger_stream_callback(logger_stream_callback, user_data)
    @ccall libuno.uno_set_logger_stream_callback(logger_stream_callback::Ptr{Cvoid},
                                                 user_data::Ptr{Cvoid})::Bool
end

function uno_reset_logger_stream()
    @ccall libuno.uno_reset_logger_stream()::Bool
end

function uno_optimize(solver, model)
    @ccall libuno.uno_optimize(solver::Ptr{Cvoid}, model::Ptr{Cvoid})::Cvoid
end

function uno_get_solver_integer_option(solver, option_name)
    @ccall libuno.uno_get_solver_integer_option(solver::Ptr{Cvoid},
                                                option_name::Cstring)::Cint
end

function uno_get_solver_double_option(solver, option_name)
    @ccall libuno.uno_get_solver_double_option(solver::Ptr{Cvoid},
                                               option_name::Cstring)::Cdouble
end

function uno_get_solver_bool_option(solver, option_name)
    @ccall libuno.uno_get_solver_bool_option(solver::Ptr{Cvoid}, option_name::Cstring)::Bool
end

function uno_get_solver_string_option(solver, option_name)
    @ccall libuno.uno_get_solver_string_option(solver::Ptr{Cvoid},
                                               option_name::Cstring)::Cstring
end

function uno_get_optimization_status(solver)
    @ccall libuno.uno_get_optimization_status(solver::Ptr{Cvoid})::Int32
end

function uno_get_solution_status(solver)
    @ccall libuno.uno_get_solution_status(solver::Ptr{Cvoid})::Int32
end

function uno_get_solution_objective(solver)
    @ccall libuno.uno_get_solution_objective(solver::Ptr{Cvoid})::Cdouble
end

function uno_get_primal_solution_component(solver, index)
    @ccall libuno.uno_get_primal_solution_component(solver::Ptr{Cvoid},
                                                    index::Int32)::Cdouble
end

function uno_get_constraint_dual_solution_component(solver, index)
    @ccall libuno.uno_get_constraint_dual_solution_component(solver::Ptr{Cvoid},
                                                             index::Int32)::Cdouble
end

function uno_get_lower_bound_dual_solution_component(solver, index)
    @ccall libuno.uno_get_lower_bound_dual_solution_component(solver::Ptr{Cvoid},
                                                              index::Int32)::Cdouble
end

function uno_get_upper_bound_dual_solution_component(solver, index)
    @ccall libuno.uno_get_upper_bound_dual_solution_component(solver::Ptr{Cvoid},
                                                              index::Int32)::Cdouble
end

function uno_get_primal_solution(solver, primal_solution)
    @ccall libuno.uno_get_primal_solution(solver::Ptr{Cvoid},
                                          primal_solution::Ptr{Cdouble})::Cvoid
end

function uno_get_constraint_dual_solution(solver, constraint_dual_solution)
    @ccall libuno.uno_get_constraint_dual_solution(solver::Ptr{Cvoid},
                                                   constraint_dual_solution::Ptr{Cdouble})::Cvoid
end

function uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)
    @ccall libuno.uno_get_lower_bound_dual_solution(solver::Ptr{Cvoid},
                                                    lower_bound_dual_solution::Ptr{Cdouble})::Cvoid
end

function uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)
    @ccall libuno.uno_get_upper_bound_dual_solution(solver::Ptr{Cvoid},
                                                    upper_bound_dual_solution::Ptr{Cdouble})::Cvoid
end

function uno_get_solution_primal_feasibility(solver)
    @ccall libuno.uno_get_solution_primal_feasibility(solver::Ptr{Cvoid})::Cdouble
end

function uno_get_solution_stationarity(solver)
    @ccall libuno.uno_get_solution_stationarity(solver::Ptr{Cvoid})::Cdouble
end

function uno_get_solution_complementarity(solver)
    @ccall libuno.uno_get_solution_complementarity(solver::Ptr{Cvoid})::Cdouble
end

function uno_get_number_iterations(solver)
    @ccall libuno.uno_get_number_iterations(solver::Ptr{Cvoid})::Int32
end

function uno_get_cpu_time(solver)
    @ccall libuno.uno_get_cpu_time(solver::Ptr{Cvoid})::Cdouble
end

function uno_get_number_objective_evaluations(solver)
    @ccall libuno.uno_get_number_objective_evaluations(solver::Ptr{Cvoid})::Int32
end

function uno_get_number_constraint_evaluations(solver)
    @ccall libuno.uno_get_number_constraint_evaluations(solver::Ptr{Cvoid})::Int32
end

function uno_get_number_objective_gradient_evaluations(solver)
    @ccall libuno.uno_get_number_objective_gradient_evaluations(solver::Ptr{Cvoid})::Int32
end

function uno_get_number_jacobian_evaluations(solver)
    @ccall libuno.uno_get_number_jacobian_evaluations(solver::Ptr{Cvoid})::Int32
end

function uno_get_number_hessian_evaluations(solver)
    @ccall libuno.uno_get_number_hessian_evaluations(solver::Ptr{Cvoid})::Int32
end

function uno_get_number_subproblem_solved_evaluations(solver)
    @ccall libuno.uno_get_number_subproblem_solved_evaluations(solver::Ptr{Cvoid})::Int32
end

function uno_destroy_model(model)
    @ccall libuno.uno_destroy_model(model::Ptr{Cvoid})::Cvoid
end

function uno_destroy_solver(solver)
    @ccall libuno.uno_destroy_solver(solver::Ptr{Cvoid})::Cvoid
end
