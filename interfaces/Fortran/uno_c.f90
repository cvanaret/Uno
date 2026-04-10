! Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
! Licensed under the MIT license. See LICENSE file in the project directory for details.

!==============================================================
! Fortran interfaces -- uno_c.f90
!==============================================================

use, intrinsic :: iso_c_binding

!---------------------------------------------
! uno_int
!---------------------------------------------
integer, parameter :: uno_int = c_int32_t

!---------------------------------------------
! Optimization sense
!---------------------------------------------
integer(uno_int), parameter :: UNO_MINIMIZE =  1
integer(uno_int), parameter :: UNO_MAXIMIZE = -1

!---------------------------------------------
! Lagrange multiplier sign convention
!---------------------------------------------
integer(uno_int), parameter :: UNO_MULTIPLIER_POSITIVE =  1
integer(uno_int), parameter :: UNO_MULTIPLIER_NEGATIVE = -1

!---------------------------------------------
! Problem type
!---------------------------------------------
character(len=2), parameter :: UNO_PROBLEM_LINEAR = "LP"
character(len=2), parameter :: UNO_PROBLEM_QUADRATIC = "QP"
character(len=3), parameter :: UNO_PROBLEM_NONLINEAR = "NLP"

!---------------------------------------------
! Base indexing style
!---------------------------------------------
integer(uno_int), parameter :: UNO_ZERO_BASED_INDEXING = 0
integer(uno_int), parameter :: UNO_ONE_BASED_INDEXING  = 1

!---------------------------------------------
! Triangular part
!---------------------------------------------
character(c_char), parameter :: UNO_LOWER_TRIANGLE = 'L'
character(c_char), parameter :: UNO_UPPER_TRIANGLE = 'U'

!---------------------------------------------
! Option type
!---------------------------------------------
integer(uno_int), parameter :: UNO_OPTION_TYPE_INTEGER   =  0
integer(uno_int), parameter :: UNO_OPTION_TYPE_DOUBLE    =  1
integer(uno_int), parameter :: UNO_OPTION_TYPE_BOOL      =  2
integer(uno_int), parameter :: UNO_OPTION_TYPE_STRING    =  3
integer(uno_int), parameter :: UNO_OPTION_TYPE_NOT_FOUND = -1

!---------------------------------------------
! Optimization status
!---------------------------------------------
integer(uno_int), parameter :: UNO_SUCCESS           = 0
integer(uno_int), parameter :: UNO_ITERATION_LIMIT   = 1
integer(uno_int), parameter :: UNO_TIME_LIMIT        = 2
integer(uno_int), parameter :: UNO_EVALUATION_ERROR  = 3
integer(uno_int), parameter :: UNO_ALGORITHMIC_ERROR = 4
integer(uno_int), parameter :: UNO_USER_TERMINATION  = 5

!---------------------------------------------
! Iterate status
!---------------------------------------------
integer(uno_int), parameter :: UNO_NOT_OPTIMAL                 = 0
integer(uno_int), parameter :: UNO_FEASIBLE_KKT_POINT          = 1
integer(uno_int), parameter :: UNO_FEASIBLE_FJ_POINT           = 2
integer(uno_int), parameter :: UNO_INFEASIBLE_STATIONARY_POINT = 3
integer(uno_int), parameter :: UNO_FEASIBLE_SMALL_STEP         = 4
integer(uno_int), parameter :: UNO_INFEASIBLE_SMALL_STEP       = 5
integer(uno_int), parameter :: UNO_UNBOUNDED                   = 6

!---------------------------------------------
! Uno version
!---------------------------------------------
integer(uno_int), parameter :: UNO_VERSION_MAJOR = 2
integer(uno_int), parameter :: UNO_VERSION_MINOR = 7
integer(uno_int), parameter :: UNO_VERSION_PATCH = 0

!---------------------------------------------
! uno_objective_callback
!---------------------------------------------
abstract interface
   function uno_objective_callback(number_variables, x, objective_value, user_data) bind(C)
      import :: uno_int, c_double, c_ptr
      integer(uno_int), value :: number_variables
      real(c_double), intent(in) :: x(*)
      real(c_double), intent(out) :: objective_value
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_objective_callback
   end function
end interface

!---------------------------------------------
! uno_constraints_callback
!---------------------------------------------
abstract interface
   function uno_constraints_callback(number_variables, number_constraints, x, constraint_values, user_data) bind(C)
      import :: uno_int, c_double, c_ptr
      integer(uno_int), value :: number_variables, number_constraints
      real(c_double), intent(in) :: x(*)
      real(c_double), intent(out) :: constraint_values(*)
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_constraints_callback
   end function
end interface

!---------------------------------------------
! uno_objective_gradient_callback
!---------------------------------------------
abstract interface
   function uno_objective_gradient_callback(number_variables, x, gradient, user_data) bind(C)
      import :: uno_int, c_double, c_ptr
      integer(uno_int), value :: number_variables
      real(c_double), intent(in) :: x(*)
      real(c_double), intent(out) :: gradient(*)
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_objective_gradient_callback
   end function
end interface

!---------------------------------------------
! uno_constraints_jacobian_callback
!---------------------------------------------
abstract interface
   function uno_constraints_jacobian_callback(number_variables, number_jacobian_nonzeros, x, jacobian_values, user_data) bind(C)
      import :: uno_int, c_double, c_ptr
      integer(uno_int), value :: number_variables, number_jacobian_nonzeros
      real(c_double), intent(in) :: x(*)
      real(c_double), intent(out) :: jacobian_values(*)
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_constraints_jacobian_callback
   end function
end interface

!---------------------------------------------
! uno_lagrangian_hessian_callback
!---------------------------------------------
abstract interface
   function uno_lagrangian_hessian_callback(number_variables, number_constraints, number_hessian_nonzeros, &
                                            x, objective_multiplier, multipliers, hessian_values, user_data) bind(C)
      import :: uno_int, c_double, c_ptr
      integer(uno_int), value :: number_variables, number_constraints, number_hessian_nonzeros
      real(c_double), intent(in) :: x(*), multipliers(*)
      real(c_double), intent(out) :: hessian_values(*)
      real(c_double), value :: objective_multiplier
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_lagrangian_hessian_callback
   end function
end interface

!---------------------------------------------
! uno_constraints_jacobian_operator_callback
!---------------------------------------------
abstract interface
   function uno_constraints_jacobian_operator_callback(number_variables, number_constraints, x, evaluate_at_x, &
                                                       vector, result, user_data) bind(C)
      import :: uno_int, c_double, c_bool, c_ptr
      integer(uno_int), value :: number_variables, number_constraints
      real(c_double), intent(in) :: x(*)
      logical(c_bool), value :: evaluate_at_x
      real(c_double), intent(in) :: vector(*)
      real(c_double), intent(out) :: result(*)
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_constraints_jacobian_operator_callback
   end function
end interface

!---------------------------------------------
! uno_constraints_jacobian_transposed_operator_callback
!---------------------------------------------
abstract interface
   function uno_constraints_jacobian_transposed_operator_callback(number_variables, number_constraints, x, &
                                                                  evaluate_at_x, vector, result, user_data) bind(C)
      import :: uno_int, c_double, c_bool, c_ptr
      integer(uno_int), value :: number_variables, number_constraints
      real(c_double), intent(in) :: x(*)
      logical(c_bool), value :: evaluate_at_x
      real(c_double), intent(in) :: vector(*)
      real(c_double), intent(out) :: result(*)
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_constraints_jacobian_transposed_operator_callback
   end function
end interface

!---------------------------------------------
! uno_lagrangian_hessian_operator_callback
!---------------------------------------------
abstract interface
   function uno_lagrangian_hessian_operator_callback(number_variables, number_constraints, x, evaluate_at_x, &
                                                     objective_multiplier, multipliers, vector, result, user_data) bind(C)
      import :: uno_int, c_double, c_bool, c_ptr
      integer(uno_int), value :: number_variables, number_constraints
      real(c_double), intent(in) :: x(*)
      logical(c_bool), value :: evaluate_at_x
      real(c_double), value :: objective_multiplier
      real(c_double), intent(in) :: multipliers(*), vector(*)
      real(c_double), intent(out) :: result(*)
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_lagrangian_hessian_operator_callback
   end function
end interface

!---------------------------------------------
! uno_notify_acceptable_iterate_callback
!---------------------------------------------
abstract interface
   subroutine uno_notify_acceptable_iterate_callback(number_variables, number_constraints, primals, &
                                                     lower_bound_multipliers, upper_bound_multipliers, &
                                                     constraint_multipliers, objective_multiplier, &
                                                     primal_feasibility_residual, stationarity_residual, &
                                                     complementarity_residual, user_data) bind(C)
      import :: uno_int, c_double, c_ptr
      integer(uno_int), value :: number_variables, number_constraints
      real(c_double), intent(in) :: primals(*)
      real(c_double), intent(in) :: lower_bound_multipliers(*)
      real(c_double), intent(in) :: upper_bound_multipliers(*)
      real(c_double), intent(in) :: constraint_multipliers(*)
      real(c_double), value :: objective_multiplier
      real(c_double), value :: primal_feasibility_residual
      real(c_double), value :: stationarity_residual
      real(c_double), value :: complementarity_residual
      type(c_ptr), value :: user_data
   end subroutine
end interface

!---------------------------------------------
! uno_termination_callback
!---------------------------------------------
abstract interface
   function uno_termination_callback(number_variables, number_constraints, primals, &
                                     lower_bound_multipliers, upper_bound_multipliers, &
                                     constraint_multipliers, objective_multiplier, &
                                     primal_feasibility_residual, stationarity_residual, &
                                     complementarity_residual, user_data) bind(C)
      import :: uno_int, c_double, c_bool, c_ptr
      integer(uno_int), value :: number_variables, number_constraints
      real(c_double), intent(in) :: primals(*)
      real(c_double), intent(in) :: lower_bound_multipliers(*)
      real(c_double), intent(in) :: upper_bound_multipliers(*)
      real(c_double), intent(in) :: constraint_multipliers(*)
      real(c_double), value :: objective_multiplier
      real(c_double), value :: primal_feasibility_residual
      real(c_double), value :: stationarity_residual
      real(c_double), value :: complementarity_residual
      type(c_ptr), value :: user_data
      logical(c_bool) :: uno_termination_callback
   end function
end interface

!---------------------------------------------
! uno_logger_stream_callback
!---------------------------------------------
abstract interface
   function uno_logger_stream_callback(buffer, length, user_data) bind(C)
      import :: uno_int, c_char, c_ptr
      character(kind=c_char), intent(in) :: buffer(*)
      integer(uno_int), value :: length
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_logger_stream_callback
   end function
end interface

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
                                jacobian_column_indices, jacobian) result(success) &
      bind(C, name="uno_set_constraints")
      import :: c_ptr, uno_int, c_double, c_bool, c_funptr
      type(c_ptr), value :: model
      integer(uno_int), value :: number_constraints, number_jacobian_nonzeros
      type(c_funptr), value :: constraint_functions, jacobian
      real(c_double) :: constraints_lower_bounds(*), constraints_upper_bounds(*)
      integer(uno_int) :: jacobian_row_indices(*), jacobian_column_indices(*)
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
                                       hessian_row_indices, hessian_column_indices, lagrangian_hessian) result(success) &
      bind(C, name="uno_set_lagrangian_hessian")
      import :: c_ptr, uno_int, c_double, c_bool, c_char, c_funptr
      type(c_ptr), value :: model
      integer(uno_int), value :: number_hessian_nonzeros
      character(c_char), value :: hessian_triangular_part
      integer(uno_int) :: hessian_row_indices(*), hessian_column_indices(*)
      type(c_funptr), value :: lagrangian_hessian
      logical(c_bool) :: success
   end function uno_set_lagrangian_hessian
end interface

!---------------------------------------------
! uno_set_lagrangian_hessian_operator
!---------------------------------------------
interface
   function uno_set_lagrangian_hessian_operator(model, lagrangian_hessian_operator) result(success) &
      bind(C, name="uno_set_lagrangian_hessian_operator")
      import :: c_ptr, c_bool, c_funptr
      type(c_ptr), value :: model
      type(c_funptr), value :: lagrangian_hessian_operator
      logical(c_bool) :: success
   end function uno_set_lagrangian_hessian_operator
end interface

!---------------------------------------------
! uno_set_lagrangian_sign_convention
!---------------------------------------------
interface
   function uno_set_lagrangian_sign_convention(model, lagrangian_sign_convention) result(success) &
      bind(C, name="uno_set_lagrangian_sign_convention")
      import :: c_ptr, uno_int, c_bool
      type(c_ptr), value :: model
      integer(uno_int), value :: lagrangian_sign_convention
      logical(c_bool) :: success
   end function uno_set_lagrangian_sign_convention
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
      real(c_double) :: initial_primal_iterate(*)
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
      real(c_double) :: initial_dual_iterate(*)
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
! uno_set_solver_callbacks
!---------------------------------------------
interface
   function uno_set_solver_callbacks(solver, notify_acceptable_iterate_callback, &
                                     termination_callback, user_data) &
                                     result(success) &
      bind(C, name="uno_set_solver_callbacks")
      import :: c_ptr, c_funptr, c_bool
      type(c_ptr), value :: solver, user_data
      type(c_funptr), value :: notify_acceptable_iterate_callback, termination_callback
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
      real(c_double) :: primal_solution(*)
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
      real(c_double) :: constraint_dual_solution(*)
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
      real(c_double) :: lower_bound_dual_solution(*)
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
      real(c_double) :: upper_bound_dual_solution(*)
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
