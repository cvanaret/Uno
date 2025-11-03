!==============================================================
!  UNO C API Fortran interfaces
!==============================================================

use, intrinsic :: iso_c_binding
implicit none

!---------------------------------------------
! uno_int
!---------------------------------------------
integer, parameter :: uno_int = c_int32_t

!---------------------------------------------
! uno_get_version
!---------------------------------------------
interface
   subroutine uno_get_version(major, minor, patch) &
      bind(C, name="uno_get_version")
      import :: uno_int
      integer(uno_int) :: major, minor, patch
   end subroutine uno_get_version
end interface

!---------------------------------------------
! uno_create_model
!---------------------------------------------
interface
   function uno_create_model(problem_type, number_variables, variables_lower_bounds, &
                             variables_upper_bounds, base_indexing) result(model) &
      bind(C, name="uno_create_model")
      import :: c_char, uno_int, c_double, c_ptr
      character(c_char), dimension(*) :: problem_type
      integer(uno_int), value :: number_variables, base_indexing
      real(c_double), dimension(*) :: variables_lower_bounds, variables_upper_bounds
      type(c_ptr) :: model
   end function uno_create_model
end interface

!---------------------------------------------
! uno_set_objective
!---------------------------------------------
interface
   function uno_set_objective(model, optimization_sense, objective_function, &
                              objective_gradient) result(success) &
      bind(C, name="uno_set_objective")
      import :: c_ptr, uno_int, c_funptr, c_bool
      type(c_ptr), value :: model
      integer(uno_int), value :: optimization_sense
      type(c_funptr), value :: objective_function
      type(c_funptr), value :: objective_gradient
      logical(c_bool) :: success
   end function uno_set_objective
end interface

!---------------------------------------------
! uno_set_constraints
!---------------------------------------------
interface
   function uno_set_constraints(model, number_constraints, constraint_functions,        &
                                constraints_lower_bounds, constraints_upper_bounds,     &
                                number_jacobian_nonzeros, jacobian_row_indices,         &
                                jacobian_column_indices, constraint_jacobian) result(success) &
      bind(C, name="uno_set_constraints")
      import :: c_ptr, uno_int, c_double, c_bool, c_funptr
      type(c_ptr), value :: model
      integer(uno_int), value :: number_constraints, number_jacobian_nonzeros
      type(c_funptr), value :: constraint_functions, constraint_jacobian
      real(c_double), dimension(*) :: constraints_lower_bounds, constraints_upper_bounds
      integer(uno_int), dimension(*) :: jacobian_row_indices, jacobian_column_indices
      logical(c_bool) :: success
   end function uno_set_constraints
end interface

!---------------------------------------------
! uno_set_jacobian_operator
!---------------------------------------------
interface
   function uno_set_jacobian_operator(model, jacobian_operator) result(success) &
      bind(C, name="uno_set_jacobian_operator")
      import :: c_ptr, c_bool, c_funptr
      type(c_ptr), value :: model
      type(c_funptr), value :: jacobian_operator
      logical(c_bool) :: success
   end function uno_set_jacobian_operator
end interface

!---------------------------------------------
! uno_set_jacobian_transposed_operator
!---------------------------------------------
interface
   function uno_set_jacobian_transposed_operator(model, jacobian_transposed_operator) result(success) &
      bind(C, name="uno_set_jacobian_transposed_operator")
      import :: c_ptr, c_bool, c_funptr
      type(c_ptr), value :: model
      type(c_funptr), value :: jacobian_transposed_operator
      logical(c_bool) :: success
   end function uno_set_jacobian_transposed_operator
end interface

!---------------------------------------------
! uno_set_lagrangian_hessian
!---------------------------------------------
interface
   function uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part, &
                                       hessian_row_indices, hessian_column_indices, lagrangian_hessian, &
                                       lagrangian_sign_convention) result(success) &
      bind(C, name="uno_set_lagrangian_hessian")
      import :: c_ptr, uno_int, c_double, c_bool, c_char, c_funptr
      type(c_ptr), value :: model
      integer(uno_int), value :: number_hessian_nonzeros
      character(c_char), value :: hessian_triangular_part
      integer(uno_int), dimension(*) :: hessian_row_indices, hessian_column_indices
      type(c_funptr), value :: lagrangian_hessian
      real(c_double), value :: lagrangian_sign_convention
      logical(c_bool) :: success
   end function uno_set_lagrangian_hessian
end interface

!---------------------------------------------
! uno_set_lagrangian_hessian_operator
!---------------------------------------------
interface
   function uno_set_lagrangian_hessian_operator(model, lagrangian_hessian_operator, &
                                                lagrangian_sign_convention) result(success) &
      bind(C, name="uno_set_lagrangian_hessian_operator")
      import :: c_ptr, uno_int, c_double, c_bool, c_funptr
      type(c_ptr), value :: model
      type(c_funptr), value :: lagrangian_hessian_operator
      real(c_double), value :: lagrangian_sign_convention
      logical(c_bool) :: success
   end function uno_set_lagrangian_hessian_operator
end interface

!---------------------------------------------
! uno_set_user_data
!---------------------------------------------
interface
   function uno_set_user_data(model, user_data) result(success) &
      bind(C, name="uno_set_user_data")
      import :: c_ptr, c_bool
      type(c_ptr), value :: model, user_data
      logical(c_bool) :: success
   end function uno_set_user_data
end interface

!---------------------------------------------
! uno_set_initial_primal_iterate_component
!---------------------------------------------
interface
   function uno_set_initial_primal_iterate_component(model, index, initial_primal_component) result(success) &
      bind(C, name="uno_set_initial_primal_iterate_component")
      import :: c_ptr, uno_int, c_double, c_bool
      type(c_ptr), value :: model
      integer(uno_int), value :: index
      real(c_double) :: initial_primal_component
      logical(c_bool) :: success
   end function uno_set_initial_primal_iterate_component
end interface

!---------------------------------------------
! uno_set_initial_dual_iterate_component
!---------------------------------------------
interface
   function uno_set_initial_dual_iterate_component(model, index, initial_dual_component) result(success) &
      bind(C, name="uno_set_initial_dual_iterate_component")
      import :: c_ptr, uno_int, c_double, c_bool
      type(c_ptr), value :: model
      integer(uno_int), value :: index
      real(c_double) :: initial_dual_component
      logical(c_bool) :: success
   end function uno_set_initial_dual_iterate_component
end interface

!---------------------------------------------
! uno_set_initial_primal_iterate
!---------------------------------------------
interface
   function uno_set_initial_primal_iterate(model, initial_primal_iterate) result(success) &
      bind(C, name="uno_set_initial_primal_iterate")
      import :: c_ptr, c_double, c_bool
      type(c_ptr), value :: model
      real(c_double), dimension(*) :: initial_primal_iterate
      logical(c_bool) :: success
   end function uno_set_initial_primal_iterate
end interface

!---------------------------------------------
! uno_set_initial_dual_iterate
!---------------------------------------------
interface
   function uno_set_initial_dual_iterate(model, initial_dual_iterate) result(success) &
      bind(C, name="uno_set_initial_dual_iterate")
      import :: c_ptr, c_double, c_bool
      type(c_ptr), value :: model
      real(c_double), dimension(*) :: initial_dual_iterate
      logical(c_bool) :: success
   end function uno_set_initial_dual_iterate
end interface

!---------------------------------------------
! uno_create_solver
!---------------------------------------------
interface
   function uno_create_solver() result(solver) &
      bind(C, name="uno_create_solver")
      import :: c_ptr
      type(c_ptr) :: solver
   end function uno_create_solver
end interface

!---------------------------------------------
! uno_set_solver_integer_option
!---------------------------------------------
interface
   function uno_set_solver_integer_option(solver, option_name, option_value) result(success) &
      bind(C, name="uno_set_solver_integer_option")
      import :: c_ptr, c_char, uno_int, c_bool
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      integer(uno_int), value :: option_value
      logical(c_bool) :: success
   end function uno_set_solver_integer_option
end interface

!---------------------------------------------
! uno_set_solver_double_option
!---------------------------------------------
interface
   function uno_set_solver_double_option(solver, option_name, option_value) result(success) &
      bind(C, name="uno_set_solver_double_option")
      import :: c_ptr, c_char, c_double, c_bool
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      real(c_double), value :: option_value
      logical(c_bool) :: success
   end function uno_set_solver_double_option
end interface

!---------------------------------------------
! uno_set_solver_bool_option
!---------------------------------------------
interface
   function uno_set_solver_bool_option(solver, option_name, option_value) result(success) &
      bind(C, name="uno_set_solver_bool_option")
      import :: c_ptr, c_char, c_bool
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      logical(c_bool), value :: option_value
      logical(c_bool) :: success
   end function uno_set_solver_bool_option
end interface

!---------------------------------------------
! uno_set_solver_string_option
!---------------------------------------------
interface
   function uno_set_solver_string_option(solver, option_name, option_value) result(success) &
      bind(C, name="uno_set_solver_string_option")
      import :: c_ptr, c_char, c_bool
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      character(c_char), dimension(*) :: option_value
      logical(c_bool) :: success
   end function uno_set_solver_string_option
end interface

!---------------------------------------------
! uno_get_solver_option_type
!---------------------------------------------
interface
   function uno_get_solver_option_type(solver, option_name) result(option_type) &
      bind(C, name="uno_get_solver_option_type")
      import :: c_ptr, c_char, uno_int
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      integer(uno_int) :: option_type
   end function uno_get_solver_option_type
end interface

!---------------------------------------------
! uno_load_solver_option_file
!---------------------------------------------
interface
   function uno_load_solver_option_file(solver, file_name) result(success) &
      bind(C, name="uno_load_solver_option_file")
      import :: c_ptr, c_char, c_bool
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: file_name
      logical(c_bool) :: success
   end function uno_load_solver_option_file
end interface

!---------------------------------------------
! uno_set_solver_preset
!---------------------------------------------
interface
   function uno_set_solver_preset(solver, preset_name) result(success) &
      bind(C, name="uno_set_solver_preset")
      import :: c_ptr, c_char, c_bool
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: preset_name
      logical(c_bool) :: success
   end function uno_set_solver_preset
end interface

!---------------------------------------------
! uno_set_solver_callbacks
!---------------------------------------------
interface
   function uno_set_solver_callbacks(solver, notify_acceptable_iterate_callback, &
                                     user_termination_callback, user_data) &
                                     result(success) &
      bind(C, name="uno_set_solver_callbacks")
      import :: c_ptr, c_funptr, c_bool
      type(c_ptr), value :: solver, user_data
      type(c_funptr), value :: notify_acceptable_iterate_callback, user_termination_callback
      logical(c_bool) :: success
   end function uno_set_solver_callbacks
end interface

!---------------------------------------------
! uno_set_logger_stream_callback
!---------------------------------------------
interface
   function uno_set_logger_stream_callback(logger_stream_callback, user_data) &
                                           result(success) &
      bind(C, name="uno_set_logger_stream_callback")
      import :: c_ptr, c_funptr, c_bool
      type(c_ptr), value :: user_data
      type(c_funptr), value :: logger_stream_callback
      logical(c_bool) :: success
   end function uno_set_logger_stream_callback
end interface

!---------------------------------------------
! uno_reset_logger_stream
!---------------------------------------------
interface
   function uno_reset_logger_stream() result(success) &
      bind(C, name="uno_reset_logger_stream")
      import :: c_bool
      logical(c_bool) :: success
   end function uno_reset_logger_stream
end interface

!---------------------------------------------
! uno_optimize
!---------------------------------------------
interface
   subroutine uno_optimize(solver, model) &
      bind(C, name="uno_optimize")
      import :: c_ptr
      type(c_ptr), value :: solver, model
   end subroutine uno_optimize
end interface

!---------------------------------------------
! uno_get_solver_integer_option
!---------------------------------------------
interface
   function uno_get_solver_integer_option(solver, option_name) result(solver_integer_option) &
      bind(C, name="uno_get_solver_integer_option")
      import :: c_ptr, c_char, uno_int
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      integer(uno_int) :: solver_integer_option
   end function uno_get_solver_integer_option
end interface

!---------------------------------------------
! uno_get_solver_double_option
!---------------------------------------------
interface
   function uno_get_solver_double_option(solver, option_name) result(solver_double_option) &
      bind(C, name="uno_get_solver_double_option")
      import :: c_ptr, c_char, c_double
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      real(c_double) :: solver_double_option
   end function uno_get_solver_double_option
end interface

!---------------------------------------------
! uno_get_solver_bool_option
!---------------------------------------------
interface
   function uno_get_solver_bool_option(solver, option_name) result(solver_bool_option) &
      bind(C, name="uno_get_solver_bool_option")
      import :: c_ptr, c_char, c_bool
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      logical(c_bool) :: solver_bool_option
   end function uno_get_solver_bool_option
end interface

!---------------------------------------------
! uno_get_solver_string_option
!---------------------------------------------
interface
   function uno_get_solver_string_option(solver, option_name) result(solver_string_option) &
      bind(C, name="uno_get_solver_string_option")
      import :: c_ptr, c_char
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      type(c_ptr) :: solver_string_option
   end function uno_get_solver_string_option
end interface

!---------------------------------------------
! uno_get_solution_status
!---------------------------------------------
interface
   function uno_get_optimization_status(solver) result(optimization_status) &
      bind(C, name="uno_get_optimization_status")
      import :: c_ptr, uno_int
      type(c_ptr), value :: solver
      integer(uno_int) :: optimization_status
   end function uno_get_optimization_status
end interface

!---------------------------------------------
! uno_get_solution_status
!---------------------------------------------
interface
   function uno_get_solution_status(solver) result(solution_status) &
      bind(C, name="uno_get_solution_status")
      import :: c_ptr, uno_int
      type(c_ptr), value :: solver
      integer(uno_int) :: solution_status
   end function uno_get_solution_status
end interface

!---------------------------------------------
! uno_get_solution_objective
!---------------------------------------------
interface
   function uno_get_solution_objective(solver) result(solution_objective) &
      bind(C, name="uno_get_solution_objective")
      import :: c_ptr, c_double
      type(c_ptr), value :: solver
      real(c_double) :: solution_objective
   end function uno_get_solution_objective
end interface

!---------------------------------------------
! uno_get_primal_solution_component
!---------------------------------------------
interface
   function uno_get_primal_solution_component(solver, index) result(primal_solution_component) &
      bind(C, name="uno_get_primal_solution_component")
      import :: c_ptr, uno_int, c_double
      type(c_ptr), value :: solver
      integer(uno_int), value :: index
      real(c_double) :: primal_solution_component
   end function uno_get_primal_solution_component
end interface

!---------------------------------------------
! uno_get_constraint_dual_solution_component
!---------------------------------------------
interface
   function uno_get_constraint_dual_solution_component(solver, index) result(constraint_dual_solution_component) &
      bind(C, name="uno_get_constraint_dual_solution_component")
      import :: c_ptr, uno_int, c_double
      type(c_ptr), value :: solver
      integer(uno_int), value :: index
      real(c_double) :: constraint_dual_solution_component
   end function uno_get_constraint_dual_solution_component
end interface

!---------------------------------------------
! uno_get_lower_bound_dual_solution_component
!---------------------------------------------
interface
   function uno_get_lower_bound_dual_solution_component(solver, index) result(lower_bound_dual_solution_component) &
      bind(C, name="uno_get_lower_bound_dual_solution_component")
      import :: c_ptr, uno_int, c_double
      type(c_ptr), value :: solver
      integer(uno_int), value :: index
      real(c_double) :: lower_bound_dual_solution_component
   end function uno_get_lower_bound_dual_solution_component
end interface

!---------------------------------------------
! uno_get_upper_bound_dual_solution_component
!---------------------------------------------
interface
   function uno_get_upper_bound_dual_solution_component(solver, index) result(upper_bound_dual_solution_component) &
      bind(C, name="uno_get_upper_bound_dual_solution_component")
      import :: c_ptr, uno_int, c_double
      type(c_ptr), value :: solver
      integer(uno_int), value :: index
      real(c_double) :: upper_bound_dual_solution_component
   end function uno_get_upper_bound_dual_solution_component
end interface

!---------------------------------------------
! uno_get_primal_solution
!---------------------------------------------
interface
   subroutine uno_get_primal_solution(solver, primal_solution) &
      bind(C, name="uno_get_primal_solution")
      import :: c_ptr, c_double
      type(c_ptr), value :: solver
      real(c_double), dimension(*) :: primal_solution
   end subroutine uno_get_primal_solution
end interface

!---------------------------------------------
! uno_get_constraint_dual_solution
!---------------------------------------------
interface
   subroutine uno_get_constraint_dual_solution(solver, constraint_dual_solution) &
      bind(C, name="uno_get_constraint_dual_solution")
      import :: c_ptr, c_double
      type(c_ptr), value :: solver
      real(c_double), dimension(*) :: constraint_dual_solution
   end subroutine uno_get_constraint_dual_solution
end interface

!---------------------------------------------
! uno_get_lower_bound_dual_solution
!---------------------------------------------
interface
   subroutine uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution) &
      bind(C, name="uno_get_lower_bound_dual_solution")
      import :: c_ptr, c_double
      type(c_ptr), value :: solver
      real(c_double), dimension(*) :: lower_bound_dual_solution
   end subroutine uno_get_lower_bound_dual_solution
end interface

!---------------------------------------------
! uno_get_upper_bound_dual_solution
!---------------------------------------------
interface
   subroutine uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution) &
      bind(C, name="uno_get_upper_bound_dual_solution")
      import :: c_ptr, c_double
      type(c_ptr), value :: solver
      real(c_double), dimension(*) :: upper_bound_dual_solution
   end subroutine uno_get_upper_bound_dual_solution
end interface

!---------------------------------------------
! uno_get_solution_primal_feasibility
!---------------------------------------------
interface
   function uno_get_solution_primal_feasibility(solver) result(solution_primal_feasibility) &
      bind(C, name="uno_get_solution_primal_feasibility")
      import :: c_ptr, c_double
      type(c_ptr), value :: solver
      real(c_double) :: solution_primal_feasibility
   end function uno_get_solution_primal_feasibility
end interface

!---------------------------------------------
! uno_get_solution_stationarity
!---------------------------------------------
interface
   function uno_get_solution_stationarity(solver) result(solution_stationarity) &
      bind(C, name="uno_get_solution_stationarity")
      import :: c_ptr, c_double
      type(c_ptr), value :: solver
      real(c_double) :: solution_stationarity
   end function uno_get_solution_stationarity
end interface

!---------------------------------------------
! uno_get_solution_complementarity
!---------------------------------------------
interface
   function uno_get_solution_complementarity(solver) result(solution_complementarity) &
      bind(C, name="uno_get_solution_complementarity")
      import :: c_ptr, c_double
      type(c_ptr), value :: solver
      real(c_double) :: solution_complementarity
   end function uno_get_solution_complementarity
end interface

!---------------------------------------------
! uno_get_number_iterations
!---------------------------------------------
interface
   function uno_get_number_iterations(solver) result(number_iterations) &
      bind(C, name="uno_get_number_iterations")
      import :: c_ptr, uno_int
      type(c_ptr), value :: solver
      integer(uno_int) :: number_iterations
   end function uno_get_number_iterations
end interface

!---------------------------------------------
! uno_get_cpu_time
!---------------------------------------------
interface
   function uno_get_cpu_time(solver) result(cpu_time) &
      bind(C, name="uno_get_cpu_time")
      import :: c_ptr, c_double
      type(c_ptr), value :: solver
      real(c_double) :: cpu_time
   end function uno_get_cpu_time
end interface

!---------------------------------------------
! uno_get_number_objective_evaluations
!---------------------------------------------
interface
   function uno_get_number_objective_evaluations(solver) result(number_objective_evaluations) &
      bind(C, name="uno_get_number_objective_evaluations")
      import :: c_ptr, uno_int
      type(c_ptr), value :: solver
      integer(uno_int) :: number_objective_evaluations
   end function uno_get_number_objective_evaluations
end interface

!---------------------------------------------
! uno_get_number_constraint_evaluations
!---------------------------------------------
interface
   function uno_get_number_constraint_evaluations(solver) result(number_constraint_evaluations) &
      bind(C, name="uno_get_number_constraint_evaluations")
      import :: c_ptr, uno_int
      type(c_ptr), value :: solver
      integer(uno_int) :: number_constraint_evaluations
   end function uno_get_number_constraint_evaluations
end interface

!---------------------------------------------
! uno_get_number_objective_gradient_evaluations
!---------------------------------------------
interface
   function uno_get_number_objective_gradient_evaluations(solver) result(number_objective_gradient_evaluations) &
      bind(C, name="uno_get_number_objective_gradient_evaluations")
      import :: c_ptr, uno_int
      type(c_ptr), value :: solver
      integer(uno_int) :: number_objective_gradient_evaluations
   end function uno_get_number_objective_gradient_evaluations
end interface

!---------------------------------------------
! uno_get_number_jacobian_evaluations
!---------------------------------------------
interface
   function uno_get_number_jacobian_evaluations(solver) result(number_jacobian_evaluations) &
      bind(C, name="uno_get_number_jacobian_evaluations")
      import :: c_ptr, uno_int
      type(c_ptr), value :: solver
      integer(uno_int) :: number_jacobian_evaluations
   end function uno_get_number_jacobian_evaluations
end interface

!---------------------------------------------
! uno_get_number_hessian_evaluations
!---------------------------------------------
interface
   function uno_get_number_hessian_evaluations(solver) result(number_hessian_evaluations) &
      bind(C, name="uno_get_number_hessian_evaluations")
      import :: c_ptr, uno_int
      type(c_ptr), value :: solver
      integer(uno_int) :: number_hessian_evaluations
   end function uno_get_number_hessian_evaluations
end interface

!---------------------------------------------
! uno_get_number_subproblem_solved_evaluations
!---------------------------------------------
interface
   function uno_get_number_subproblem_solved_evaluations(solver) result(number_subproblem_solved_evaluations) &
      bind(C, name="uno_get_number_subproblem_solved_evaluations")
      import :: c_ptr, uno_int
      type(c_ptr), value :: solver
      integer(uno_int) :: number_subproblem_solved_evaluations
   end function uno_get_number_subproblem_solved_evaluations
end interface

!---------------------------------------------
! uno_destroy_model
!---------------------------------------------
interface
   subroutine uno_destroy_model(model) &
      bind(C, name="uno_destroy_model")
      import :: c_ptr
      type(c_ptr), value :: model
   end subroutine uno_destroy_model
end interface

!---------------------------------------------
! uno_destroy_solver
!---------------------------------------------
interface
   subroutine uno_destroy_solver(solver) &
      bind(C, name="uno_destroy_solver")
      import :: c_ptr
      type(c_ptr), value :: solver
   end subroutine uno_destroy_solver
end interface
