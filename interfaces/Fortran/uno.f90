!==============================================================
!  UNO C API Fortran interfaces
!==============================================================

use, intrinsic :: iso_c_binding
implicit none

!---------------------------------------------
! uno_get_version
!---------------------------------------------
interface
   subroutine uno_get_version(major, minor, patch) bind(C, name="uno_get_version")
      import :: c_int
      integer(c_int) :: major, minor, patch
   end subroutine uno_get_version
end interface

!---------------------------------------------
! uno_create_model
!---------------------------------------------
interface
   function uno_create_model(problem_type, number_variables, variables_lower_bounds, variables_upper_bounds, base_indexing) result(model) bind(C, name="uno_create_model")
      import :: c_char, c_int, c_double, c_ptr
      character(c_char), value :: problem_type
      integer(c_int), value :: number_variables, base_indexing
      real(c_double), dimension(*) :: variables_lower_bounds, variables_upper_bounds
      type(c_ptr) :: model
   end function uno_create_model
end interface

!---------------------------------------------
! uno_set_objective
!---------------------------------------------
interface
   function uno_set_objective(model, optimization_sense, objective_function, objective_gradient) result(success) bind(C, name="uno_set_objective")
      import :: c_ptr, c_int, c_funptr, c_bool
      type(c_ptr), value :: model
      integer(c_int), value :: optimization_sense
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
      import :: c_ptr, c_int, c_double, c_bool, c_funptr
      type(c_ptr), value :: model
      integer(c_int), value :: number_constraints, number_jacobian_nonzeros
      type(c_funptr), value :: constraint_functions, constraint_jacobian
      real(c_double), dimension(*) :: constraints_lower_bounds, constraints_upper_bounds
      integer(c_int), dimension(*) :: jacobian_row_indices, jacobian_column_indices
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
      import :: c_ptr, c_int, c_double, c_bool, c_char, c_funptr
      type(c_ptr), value :: model
      integer(c_int), value :: number_hessian_nonzeros
      character(c_char), value :: hessian_triangular_part
      integer(c_int), dimension(*) :: hessian_row_indices, hessian_column_indices
      type(c_funptr), value :: lagrangian_hessian
      real(c_double), value :: lagrangian_sign_convention
      logical(c_bool) :: success
   end function uno_set_lagrangian_hessian
end interface

!---------------------------------------------
! uno_set_lagrangian_hessian_operator
!---------------------------------------------
interface
   function uno_set_lagrangian_hessian_operator(model, number_hessian_nonzeros, lagrangian_hessian_operator, &
                                                lagrangian_sign_convention) result(success) &
      bind(C, name="uno_set_lagrangian_hessian_operator")
      import :: c_ptr, c_int, c_double, c_bool, c_funptr
      type(c_ptr), value :: model
      integer(c_int), value :: number_hessian_nonzeros
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
      import :: c_ptr, c_int, c_double, c_bool
      type(c_ptr), value :: model
      integer(c_int), value :: index
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
      import :: c_ptr, c_int, c_double, c_bool
      type(c_ptr), value :: model
      integer(c_int), value :: index
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
   function uno_create_solver() result(solver) bind(C, name="uno_create_solver")
      import :: c_ptr
      type(c_ptr) :: solver
   end function uno_create_solver
end interface

!---------------------------------------------
! uno_set_solver_integer_option
!---------------------------------------------
interface
   subroutine uno_set_solver_integer_option(solver, option_name, option_value) bind(C, name="uno_set_solver_integer_option")
      import :: c_ptr, c_char, c_int
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      integer(c_int), value :: option_value
   end subroutine uno_set_solver_integer_option
end interface

!---------------------------------------------
! uno_set_solver_double_option
!---------------------------------------------
interface
   subroutine uno_set_solver_double_option(solver, option_name, option_value) bind(C, name="uno_set_solver_double_option")
      import :: c_ptr, c_char, c_double
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      real(c_double), value :: option_value
   end subroutine uno_set_solver_double_option
end interface

!---------------------------------------------
! uno_set_solver_bool_option
!---------------------------------------------
interface
   subroutine uno_set_solver_bool_option(solver, option_name, option_value) bind(C, name="uno_set_solver_bool_option")
      import :: c_ptr, c_char, c_bool
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      logical(c_bool), value :: option_value
   end subroutine uno_set_solver_bool_option
end interface

!---------------------------------------------
! uno_set_solver_string_option
!---------------------------------------------
interface
   subroutine uno_set_solver_string_option(solver, option_name, option_value) bind(C, name="uno_set_solver_string_option")
      import :: c_ptr, c_char
      type(c_ptr), value :: solver
      character(c_char), dimension(*) :: option_name
      character(c_char), dimension(*) :: option_value
   end subroutine uno_set_solver_string_option
end interface

! // gets the type of a particular option in the Uno solver.
! // takes as input the name of the option.
! // the possible types are integer, double, bool and string.
! int32_t uno_get_solver_option_type(void* solver, const char* option_name);

! // [optional] loads the options from a given option file.
! // takes as input the name of the option file.
! void uno_load_solver_option_file(void* solver, const char* file_name);

! // sets a particular preset in the Uno solver.
! void uno_set_solver_preset(void* solver, const char* preset_name);

! // [optional]
! // sets the user callbacks for solver.
! void uno_set_solver_callbacks(void* solver, NotifyAcceptableIterateUserCallback notify_acceptable_iterate_callback,
!    TerminationUserCallback user_termination_callback, void* user_data);

! // [optional]
! // sets the logger stream callback.
! void uno_set_logger_stream_callback(LoggerStreamUserCallback logger_stream_callback, void* user_data);

! // [optional]
! // resets the logger stream to the standard output
! void uno_reset_logger_stream();

!---------------------------------------------
! uno_optimize
!---------------------------------------------
interface
   subroutine uno_optimize(solver, model) bind(C, name="uno_optimize")
      import :: c_ptr
      type(c_ptr), value :: solver, model
   end subroutine uno_optimize
end interface

! // gets the value of a given double option.
! // takes as inputs the name of the option.
! // the possible types are integer, unsigned integer, double, bool and string.
! int uno_get_solver_integer_option(void* solver, const char* option_name);
! size_t uno_get_solver_unsigned_integer_option(void* solver, const char* option_name);
! double uno_get_solver_double_option(void* solver, const char* option_name);
! bool uno_get_solver_bool_option(void* solver, const char* option_name);
! const char* uno_get_solver_string_option(void* solver, const char* option_name);

! // gets the optimization status (once the model was solved)
! int32_t uno_get_optimization_status(void* solver);

! // gets the iterate status (once the model was solved)
! int32_t uno_get_solution_status(void* solver);

! // gets the objective value at the solution (once the model was solved)
! double uno_get_solution_objective(void* solver);

! // gets one component of the primal solution (once the model was solved)
! double uno_get_primal_solution_component(void* solver, int32_t index);

! // gets one component of the constraint dual solution (once the model was solved)
! double uno_get_constraint_dual_solution_component(void* solver, int32_t index);

! // gets one component of the lower bound dual solution (once the model was solved)
! double uno_get_lower_bound_dual_solution_component(void* solver, int32_t index);

! // gets one component of the upper bound dual solution (once the model was solved)
! double uno_get_upper_bound_dual_solution_component(void* solver, int32_t index);

! // gets the primal solution (once the model was solved)
! void uno_get_primal_solution(void* solver, double* primal_solution);

! // gets the dual solution associated with the constraints (once the model was solved)
! void uno_get_constraint_dual_solution(void* solver, double* constraint_dual_solution);

! // gets the dual solution associated with the lower bounds (once the model was solved)
! void uno_get_lower_bound_dual_solution(void* solver, double* lower_bound_dual_solution);

! // gets the dual solution associated with the upper bounds (once the model was solved)
! void uno_get_upper_bound_dual_solution(void* solver, double* upper_bound_dual_solution);

! // gets the primal feasibility residual at the solution (once the model was solved)
! double uno_get_solution_primal_feasibility(void* solver);

! // gets the stationarity residual at the solution (once the model was solved)
! double uno_get_solution_stationarity(void* solver);

! // gets the complementarity residual at the solution (once the model was solved)
! double uno_get_solution_complementarity(void* solver);

! // gets the number of outer iterations required by the solver (once the model was solved)
! int32_t uno_get_number_iterations(void* solver);

! // gets the CPU time required by the solver (once the model was solved)
! double uno_get_cpu_time(void* solver);

! // gets the number of objective evaluations required by the solver (once the model was solved)
! int32_t uno_get_number_objective_evaluations(void* solver);

! // gets the number of constraint evaluations required by the solver (once the model was solved)
! int32_t uno_get_number_constraint_evaluations(void* solver);

! // gets the number of objective gradient evaluations required by the solver (once the model was solved)
! int32_t uno_get_number_objective_gradient_evaluations(void* solver);

! // gets the number of constraint Jacobian evaluations required by the solver (once the model was solved)
! int32_t uno_get_number_jacobian_evaluations(void* solver);

! // gets the number of Lagrangian Hessian evaluations required by the solver (once the model was solved)
! int32_t uno_get_number_hessian_evaluations(void* solver);

! // gets the number of subproblems solved by the solver (once the model was solved)
! int32_t uno_get_number_subproblem_solved_evaluations(void* solver);

!---------------------------------------------
! uno_destroy_model
!---------------------------------------------
interface
   subroutine uno_destroy_model(model) bind(C, name="uno_destroy_model")
      import :: c_ptr
      type(c_ptr), value :: model
   end subroutine uno_destroy_model
end interface

!---------------------------------------------
! uno_destroy_solver
!---------------------------------------------
interface
   subroutine uno_destroy_solver(solver) bind(C, name="uno_destroy_solver")
      import :: c_ptr
      type(c_ptr), value :: solver
   end subroutine uno_destroy_solver
end interface
