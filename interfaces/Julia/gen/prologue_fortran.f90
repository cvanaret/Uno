use, intrinsic :: iso_c_binding

!---------------------------------------------
! uno_int
!---------------------------------------------
!--- UNO_INT_KIND ---
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
character(len=2), parameter :: UNO_PROBLEM_LINEAR    = "LP"
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
integer(uno_int), parameter :: UNO_VERSION_PATCH = 2

!---------------------------------------------
! uno_objective_callback
!---------------------------------------------
abstract interface
   function uno_objective_callback(number_variables, x, objective_value, user_data) &
      bind(C)
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
   function uno_constraints_callback(number_variables, number_constraints, x, &
                                     constraint_values, user_data) &
      bind(C)
      import :: uno_int, c_double, c_ptr
      integer(uno_int), value :: number_variables
      integer(uno_int), value :: number_constraints
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
   function uno_objective_gradient_callback(number_variables, x, gradient, user_data) &
      bind(C)
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
   function uno_constraints_jacobian_callback(number_variables, &
                                              number_jacobian_nonzeros, x, &
                                              jacobian_values, user_data) &
      bind(C)
      import :: uno_int, c_double, c_ptr
      integer(uno_int), value :: number_variables
      integer(uno_int), value :: number_jacobian_nonzeros
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
   function uno_lagrangian_hessian_callback(number_variables, number_constraints, &
                                            number_hessian_nonzeros, x, &
                                            objective_multiplier, multipliers, &
                                            hessian_values, user_data) &
      bind(C)
      import :: uno_int, c_double, c_ptr
      integer(uno_int), value :: number_variables
      integer(uno_int), value :: number_constraints
      integer(uno_int), value :: number_hessian_nonzeros
      real(c_double), intent(in) :: x(*)
      real(c_double), value :: objective_multiplier
      real(c_double), intent(in) :: multipliers(*)
      real(c_double), intent(out) :: hessian_values(*)
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_lagrangian_hessian_callback
   end function
end interface

!---------------------------------------------
! uno_constraints_jacobian_operator_callback
!---------------------------------------------
abstract interface
   function uno_constraints_jacobian_operator_callback(number_variables, &
                                                       number_constraints, x, &
                                                       evaluate_at_x, vector, &
                                                       result, user_data) &
      bind(C)
      import :: uno_int, c_double, c_bool, c_ptr
      integer(uno_int), value :: number_variables
      integer(uno_int), value :: number_constraints
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
   function uno_constraints_jacobian_transposed_operator_callback(number_variables, &
                                                                  number_constraints, x, &
                                                                  evaluate_at_x, vector, &
                                                                  result, user_data) &
      bind(C)
      import :: uno_int, c_double, c_bool, c_ptr
      integer(uno_int), value :: number_variables
      integer(uno_int), value :: number_constraints
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
   function uno_lagrangian_hessian_operator_callback(number_variables, &
                                                     number_constraints, x, &
                                                     evaluate_at_x, &
                                                     objective_multiplier, &
                                                     multipliers, vector, result, &
                                                     user_data) &
      bind(C)
      import :: uno_int, c_double, c_bool, c_ptr
      integer(uno_int), value :: number_variables
      integer(uno_int), value :: number_constraints
      real(c_double), intent(in) :: x(*)
      logical(c_bool), value :: evaluate_at_x
      real(c_double), value :: objective_multiplier
      real(c_double), intent(in) :: multipliers(*)
      real(c_double), intent(in) :: vector(*)
      real(c_double), intent(out) :: result(*)
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_lagrangian_hessian_operator_callback
   end function
end interface

!---------------------------------------------
! uno_notify_acceptable_iterate_callback
!---------------------------------------------
abstract interface
   subroutine uno_notify_acceptable_iterate_callback(number_variables, &
                                                     number_constraints, primals, &
                                                     lower_bound_multipliers, &
                                                     upper_bound_multipliers, &
                                                     constraint_multipliers, &
                                                     objective_multiplier, &
                                                     primal_feasibility_residual, &
                                                     stationarity_residual, &
                                                     complementarity_residual, &
                                                     user_data) &
      bind(C)
      import :: uno_int, c_double, c_ptr
      integer(uno_int), value :: number_variables
      integer(uno_int), value :: number_constraints
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
                                     complementarity_residual, user_data) &
      bind(C)
      import :: uno_int, c_double, c_ptr
      integer(uno_int), value :: number_variables
      integer(uno_int), value :: number_constraints
      real(c_double), intent(in) :: primals(*)
      real(c_double), intent(in) :: lower_bound_multipliers(*)
      real(c_double), intent(in) :: upper_bound_multipliers(*)
      real(c_double), intent(in) :: constraint_multipliers(*)
      real(c_double), value :: objective_multiplier
      real(c_double), value :: primal_feasibility_residual
      real(c_double), value :: stationarity_residual
      real(c_double), value :: complementarity_residual
      type(c_ptr), value :: user_data
      logical(uno_int) :: uno_termination_callback
   end function
end interface

!---------------------------------------------
! uno_logger_stream_callback
!---------------------------------------------
abstract interface
   function uno_logger_stream_callback(buffer, length, user_data) &
      bind(C)
      import :: c_char, uno_int, c_ptr
      character(c_char), intent(in) :: buffer(*)
      integer(uno_int), value :: length
      type(c_ptr), value :: user_data
      integer(uno_int) :: uno_logger_stream_callback
   end function
end interface
