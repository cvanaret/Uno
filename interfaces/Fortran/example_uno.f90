program example_uno
    include 'uno_c.f90'
    integer(uno_int) :: major, minor, patch
    integer(uno_int), parameter :: number_variables = 2, number_constraints = 2
    integer(uno_int), parameter :: number_jacobian_nonzeros = 4, number_hessian_nonzeros = 3
    real(c_double), dimension(number_variables) :: x0, variables_lower_bounds, variables_upper_bounds
    real(c_double), dimension(number_constraints) :: y0, constraints_lower_bounds, constraints_upper_bounds
    integer(uno_int), dimension(number_jacobian_nonzeros) :: jacobian_row_indices, jacobian_column_indices
    integer(uno_int), dimension(number_hessian_nonzeros) :: hessian_row_indices, hessian_column_indices
    type(c_ptr) :: model, solver
    real(c_double) :: solution_objective
    real(c_double), dimension(number_variables) :: primal_solution, lower_bound_dual_solution, upper_bound_dual_solution
    real(c_double), dimension(number_constraints) :: constraint_dual_solution, solution_constraints
    real(c_double) :: solution_primal_feasibility, solution_stationarity, solution_complementarity
    integer(uno_int) :: optimization_status, iterate_status, option_type, status
    logical(c_bool) :: success, print_solution = .true.
    integer(uno_int) :: number_iterations, max_iterations = 1000
    real(c_double) ::  cpu_time, dual_tolerance, primal_tolerance = 1.0d-6
    character(:), allocatable :: logger, option_file
    character(:), allocatable :: method_description
    character(len=32) :: preset = "filtersqp"  ! can be overridden by the first command-line argument
    character(len=*), parameter :: hessian_model = "exact"
    character(len=*), parameter :: problem_type = UNO_PROBLEM_NONLINEAR
    integer(uno_int), parameter :: base_indexing = UNO_ONE_BASED_INDEXING
    integer(uno_int), parameter :: optimization_sense = UNO_MINIMIZE
    character(len=1), parameter :: hessian_triangular_part = UNO_LOWER_TRIANGLE
    integer(uno_int), parameter :: lagrangian_sign_convention = UNO_MULTIPLIER_NEGATIVE
    type(c_funptr) :: objective, gradient, constraints, jacobian, lagrangian_hessian
    type(c_funptr) :: jacobian_operator, jacobian_transposed_operator, lagrangian_hessian_operator
    procedure(uno_objective_callback), pointer :: objective_callback => null()
    procedure(uno_objective_gradient_callback), pointer :: objective_gradient_callback => null()
    procedure(uno_constraints_callback), pointer :: constraints_callback => null()
    procedure(uno_constraints_jacobian_callback), pointer :: constraints_jacobian_callback => null()
    procedure(uno_lagrangian_hessian_callback), pointer :: lagrangian_hessian_callback => null()
    procedure(uno_constraints_jacobian_operator_callback), pointer :: constraints_jacobian_operator_callback => null()
    procedure(uno_constraints_jacobian_transposed_operator_callback), &
        pointer :: constraints_jacobian_transposed_operator_callback => null()
    procedure(uno_lagrangian_hessian_operator_callback), pointer :: lagrangian_hessian_operator_callback => null()

    !---------------------------------------------------
    ! Versions
    !---------------------------------------------------
    call uno_get_version(major, minor, patch)
    print *, 'Uno version: ', major, '.', minor, '.', patch

    !---------------------------------------------------
    ! Variables
    !---------------------------------------------------
    variables_lower_bounds = [-HUGE(1.0d0), -HUGE(1.0d0)]
    variables_upper_bounds = [0.5d0, HUGE(1.0d0)]
    x0 = [-2.0d0, 1.0d0]

    !---------------------------------------------------
    ! Constraints
    !---------------------------------------------------
    constraints_lower_bounds = [1.0d0, 0.0d0]
    constraints_upper_bounds = [HUGE(1.0d0), HUGE(1.0d0)]
    y0 = [0.0d0, 0.0d0]

    !---------------------------------------------------
    ! Jacobian
    !---------------------------------------------------
    jacobian_row_indices = [1, 2, 1, 2]
    jacobian_column_indices = [1, 1, 2, 2]

    !---------------------------------------------------
    ! Hessian
    !---------------------------------------------------
    hessian_row_indices = [1, 2, 2]
    hessian_column_indices = [1, 1, 2]

    !---------------------------------------------------
    ! Callbacks for Uno
    !---------------------------------------------------
    objective_callback => objective_hs15
    objective_gradient_callback  => gradient_hs15
    constraints_callback => constraints_hs15
    constraints_jacobian_callback => jacobian_hs15
    constraints_jacobian_operator_callback => jacobian_operator_hs15
    constraints_jacobian_transposed_operator_callback => jacobian_transposed_operator_hs15
    lagrangian_hessian_callback => lagrangian_hessian_hs15
    lagrangian_hessian_operator_callback => lagrangian_hessian_operator_hs15

    objective = c_funloc(objective_callback)
    gradient = c_funloc(objective_gradient_callback)
    constraints = c_funloc(constraints_callback)
    jacobian = c_funloc(constraints_jacobian_callback)
    jacobian_operator = c_funloc(constraints_jacobian_operator_callback)
    jacobian_transposed_operator = c_funloc(constraints_jacobian_transposed_operator_callback)
    lagrangian_hessian = c_funloc(lagrangian_hessian_callback)
    lagrangian_hessian_operator = c_funloc(lagrangian_hessian_operator_callback)

    !---------------------------------------------------
    ! Model creation
    !---------------------------------------------------
    model = uno_create_model(problem_type, number_variables, variables_lower_bounds, variables_upper_bounds, base_indexing)

    success = uno_set_model_name(model, "hs015")
    call check_success(success, "uno_set_model_name")
    success = uno_set_objective(model, optimization_sense, objective, gradient)
    call check_success(success, "uno_set_objective")
    success = uno_set_constraints(model, number_constraints, constraints, constraints_lower_bounds, constraints_upper_bounds, &
                                  number_jacobian_nonzeros, jacobian_row_indices, jacobian_column_indices, jacobian)
    call check_success(success, "uno_set_constraints")
    success = uno_set_jacobian_operator(model, jacobian_operator)
    call check_success(success, "uno_set_jacobian_operator")
    success = uno_set_jacobian_transposed_operator(model, jacobian_transposed_operator)
    call check_success(success, "uno_set_jacobian_transposed_operator")
    ! success = uno_set_lagrangian_hessian_operator(model, lagrangian_hessian_operator)
    ! success = uno_set_lagrangian_sign_convention(model, lagrangian_sign_convention)
    success = uno_set_initial_primal_iterate(model, x0)
    call check_success(success, "uno_set_initial_primal_iterate")
    success = uno_set_initial_dual_iterate(model, y0)
    call check_success(success, "uno_set_initial_dual_iterate")

    !---------------------------------------------------
    ! Solver creation
    !---------------------------------------------------
    solver = uno_create_solver()
    status = uno_set_solver_integer_option(solver, "max_iterations", max_iterations)
    call check_status(status, UNO_OPTIONS_SUCCESS, "uno_set_solver_integer_option")
    status = uno_set_solver_double_option(solver, "primal_tolerance", primal_tolerance)
    call check_status(status, UNO_OPTIONS_SUCCESS, "uno_set_solver_double_option")
    status = uno_set_solver_bool_option(solver, "print_solution", print_solution)
    call check_status(status, UNO_OPTIONS_SUCCESS, "uno_set_solver_bool_option")
    status = uno_set_solver_string_option(solver, "hessian_model", hessian_model)
    call check_status(status, UNO_OPTIONS_SUCCESS, "uno_set_solver_string_option")
    if (command_argument_count() >= 1) then
        call get_command_argument(1, preset)
    end if
    print *, 'preset = ', trim(preset)
    status = uno_set_solver_preset(solver, trim(preset))
    call check_status(status, UNO_OPTIONS_SUCCESS, "uno_set_solver_preset")

    ! invalid options are reported with status codes and leave the options unchanged
    status = uno_set_solver_string_option(solver, "logger", "SOLENT")
    call check_status(status, UNO_OPTIONS_INVALID_VALUE, "uno_set_solver_string_option (invalid value)")
    status = uno_set_solver_integer_option(solver, "max_iteration", max_iterations)
    call check_status(status, UNO_OPTIONS_NOT_FOUND, "uno_set_solver_integer_option (unknown option)")
    status = uno_set_solver_double_option(solver, "max_iterations", primal_tolerance)
    call check_status(status, UNO_OPTIONS_WRONG_TYPE, "uno_set_solver_double_option (wrong type)")
    status = uno_set_solver_preset(solver, "snopt")
    call check_status(status, UNO_OPTIONS_INVALID_VALUE, "uno_set_solver_preset (invalid preset)")
    status = uno_load_solver_option_file(solver, "this_option_file_does_not_exist.opt")
    call check_status(status, UNO_OPTIONS_INVALID_FILE, "uno_load_solver_option_file (missing file)")

    ! no result before the solve
    status = uno_get_optimization_status(solver, optimization_status)
    call check_status(status, UNO_RESULT_NOT_AVAILABLE, "uno_get_optimization_status (before the solve)")

    !---------------------------------------------------
    ! Run 1: solve with no Hessian. Uno defaults to L-BFGS Hessian for NLPs
    !---------------------------------------------------
    status = uno_optimize(solver, model)
    call check_status(status, UNO_OPTIMIZE_SUCCESS, "uno_optimize")
    status = uno_get_solution_objective(solver, solution_objective)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_solution_objective")
    print *, 'Solution objective = ', solution_objective

    !---------------------------------------------------
    ! Run 2: solve with exact Hessian
    !---------------------------------------------------
    success = uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part, hessian_row_indices, &
                                         hessian_column_indices, lagrangian_hessian)
    call check_success(success, "uno_set_lagrangian_hessian")
    success = uno_set_lagrangian_sign_convention(model, lagrangian_sign_convention)
    call check_success(success, "uno_set_lagrangian_sign_convention")
    status = uno_optimize(solver, model)
    call check_status(status, UNO_OPTIMIZE_SUCCESS, "uno_optimize")

    !---------------------------------------------------
    ! Get method description
    !---------------------------------------------------
    status = uno_get_method_description(solver, method_description)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_method_description")
    print *, 'method description = ', method_description
    ! the second run must use the exact Hessian set above (not fall back to L-BFGS)
    if (index(method_description, 'exact Hessian') == 0) then
        error stop 'The exact Hessian was not used in run 2'
    end if

    !---------------------------------------------------
    ! Get options
    !---------------------------------------------------
    status = uno_get_solver_integer_option(solver, "max_iterations", max_iterations)
    call check_status(status, UNO_OPTIONS_SUCCESS, "uno_get_solver_integer_option")
    print *, 'max_iterations = ', max_iterations
    if (max_iterations /= 1000) error stop 'max_iterations was modified by a rejected option'

    status = uno_get_solver_double_option(solver, "dual_tolerance", dual_tolerance)
    call check_status(status, UNO_OPTIONS_SUCCESS, "uno_get_solver_double_option")
    print *, 'dual_tolerance = ', dual_tolerance

    status = uno_get_solver_bool_option(solver, "print_solution", print_solution)
    call check_status(status, UNO_OPTIONS_SUCCESS, "uno_get_solver_bool_option")
    print *, 'print_solution = ', print_solution

    status = uno_get_solver_string_option(solver, "logger", logger)
    call check_status(status, UNO_OPTIONS_SUCCESS, "uno_get_solver_string_option")
    print *, 'logger = ', logger
    if (logger /= "INFO") error stop 'logger was modified by a rejected option'

    ! errors of the getters
    status = uno_get_solver_integer_option(solver, "max_iteration", max_iterations)
    call check_status(status, UNO_OPTIONS_NOT_FOUND, "uno_get_solver_integer_option (unknown option)")
    status = uno_get_solver_double_option(solver, "max_iterations", dual_tolerance)
    call check_status(status, UNO_OPTIONS_WRONG_TYPE, "uno_get_solver_double_option (wrong type)")
    status = uno_get_solver_string_option(solver, "option_file", option_file)
    call check_status(status, UNO_OPTIONS_NOT_SET, "uno_get_solver_string_option (option without value)")
    if (allocated(option_file)) error stop 'option_file was written although it has no value'

    option_type = uno_get_solver_option_type(solver, "linear_solver")
    print *, 'option_type for the option linear_solver = ', option_type

    !---------------------------------------------------
    ! Get solution
    !---------------------------------------------------
    status = uno_get_number_iterations(solver, number_iterations)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_number_iterations")
    print *, 'number_iterations = ', number_iterations

    status = uno_get_cpu_time(solver, cpu_time)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_cpu_time")
    print *, 'cpu_time = ', cpu_time

    status = uno_get_optimization_status(solver, optimization_status)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_optimization_status")
    print *, 'optimization_status = ', optimization_status

    status = uno_get_solution_status(solver, iterate_status)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_solution_status")
    print *, 'iterate_status = ', iterate_status

    status = uno_get_solution_objective(solver, solution_objective)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_solution_objective")
    print *, 'Solution objective = ', solution_objective

    status = uno_get_primal_solution(solver, primal_solution)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_primal_solution")
    print *, 'Primal solution = ', primal_solution

    status = uno_get_constraint_dual_solution(solver, constraint_dual_solution)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_constraint_dual_solution")
    print *, 'Constraint dual solution = ', constraint_dual_solution

    status = uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_lower_bound_dual_solution")
    print *, 'Lower bound dual solution = ', lower_bound_dual_solution

    status = uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_upper_bound_dual_solution")
    print *, 'Upper bound dual solution = ', upper_bound_dual_solution
    
    status = uno_get_solution_constraints(solver, solution_constraints)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_solution_constraints")
    print *, 'Solution constraints = ', solution_constraints

    status = uno_get_solution_primal_feasibility(solver, solution_primal_feasibility)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_solution_primal_feasibility")
    print *, 'Primal feasibility at solution = ', solution_primal_feasibility

    status = uno_get_solution_stationarity(solver, solution_stationarity)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_solution_stationarity")
    print *, 'Stationarity at solution = ', solution_stationarity

    status = uno_get_solution_complementarity(solver, solution_complementarity)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_solution_complementarity")

    ! components of the solution (1-based indexing)
    status = uno_get_primal_solution_component(solver, 2, solution_objective)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_primal_solution_component")
    if (solution_objective /= primal_solution(2)) error stop 'Wrong primal solution component'
    status = uno_get_primal_solution_component(solver, 3, solution_objective)
    call check_status(status, UNO_RESULT_INVALID_INDEX, "uno_get_primal_solution_component (invalid index)")
    status = uno_get_solution_objective(solver, solution_objective)
    call check_status(status, UNO_RESULT_SUCCESS, "uno_get_solution_objective")
    print *, 'Complementarity at solution = ', solution_complementarity

    !---------------------------------------------------
    ! Check the solution (hs015: f* = 306.5 at x* = (0.5, 2))
    !---------------------------------------------------
    if (optimization_status /= UNO_SUCCESS) error stop 'The optimization did not succeed'
    ! relative tolerance on the objective: an interior-point method stops close to, not at, the optimum
    if (abs(solution_objective - 306.5d0) > 1.0d-6 * 306.5d0) error stop 'Wrong solution objective'
    if (any(abs(primal_solution - [0.5d0, 2.0d0]) > 1.0d-6)) error stop 'Wrong primal solution'

    !---------------------------------------------------
    ! Cleanup
    !---------------------------------------------------
    deallocate(logger)
    deallocate(method_description)
    call uno_destroy_solver(solver)
    call uno_destroy_model(model)

contains
    include 'uno_fortran.f90'

    ! Stop if a call to the Uno API was rejected
    subroutine check_success(success, function_name)
        logical(c_bool), intent(in) :: success
        character(len=*), intent(in) :: function_name
        if (.not. success) then
            print *, 'Error: ', function_name, ' failed'
            error stop 1
        end if
    end subroutine check_success

    subroutine check_status(status, expected_status, function_name)
        integer(uno_int), intent(in) :: status, expected_status
        character(len=*), intent(in) :: function_name
        if (status /= expected_status) then
            print *, 'Error: ', function_name, ' returned ', status, ' instead of ', expected_status
            error stop 1
        end if
    end subroutine check_status

    ! Objective
    function objective_hs15(number_variables, x, objective_value, user_data) result(res) &
        bind(C)
        integer(uno_int), value :: number_variables
        real(c_double), intent(in) :: x(*)
        real(c_double), intent(out) :: objective_value
        type(c_ptr), value :: user_data
        integer(uno_int) :: res

        objective_value = 100.0d0 * (x(2) - x(1)**2)**2 + (1.0d0 - x(1))**2
        res = 0
    end function objective_hs15

    ! Gradient
    function gradient_hs15(number_variables, x, gradient, user_data) result(res) &
        bind(C)
        integer(uno_int), value :: number_variables
        real(c_double), intent(in) :: x(*)
        real(c_double), intent(out) :: gradient(*)
        type(c_ptr), value :: user_data
        integer(uno_int) :: res

        gradient(1) = 400.0d0 * x(1)**3 - 400.0d0 * x(1) * x(2) + 2.0d0 * x(1) - 2.0d0
        gradient(2) = 200.0d0 * (x(2) - x(1)**2)
        res = 0
    end function gradient_hs15

    ! Constraints
    function constraints_hs15(number_variables, number_constraints, x, constraint_values, user_data) result(res) &
        bind(C)
        integer(uno_int), value :: number_variables, number_constraints
        real(c_double), intent(in) :: x(*)
        real(c_double), intent(out) :: constraint_values(*)
        type(c_ptr), value :: user_data
        integer(uno_int) :: res

        constraint_values(1) = x(1) * x(2)
        constraint_values(2) = x(1) + x(2)**2
        res = 0
    end function constraints_hs15

    ! Jacobian
    function jacobian_hs15(number_variables, number_jacobian_nonzeros, x, jacobian_values, user_data) result(res) &
        bind(C)
        integer(uno_int), value :: number_variables, number_jacobian_nonzeros
        real(c_double), intent(in) :: x(*)
        real(c_double), intent(out) :: jacobian_values(*)
        type(c_ptr), value :: user_data
        integer(uno_int) :: res

        jacobian_values(1) = x(2)
        jacobian_values(2) = 1.0d0
        jacobian_values(3) = x(1)
        jacobian_values(4) = 2.0d0 * x(2)
        res = 0
    end function jacobian_hs15

    ! Jacobian operator
    function jacobian_operator_hs15(number_variables, number_constraints, x, evaluate_at_x, &
                                    vector, result, user_data) result(res) &
        bind(C)
        integer(uno_int), value :: number_variables, number_constraints
        real(c_double), intent(in) :: x(*)
        logical(c_bool), value :: evaluate_at_x
        real(c_double), intent(in) :: vector(*)
        real(c_double), intent(out) :: result(*)
        type(c_ptr), value :: user_data
        integer(uno_int) :: res

        result(1) = x(2) * vector(1) + 1.0d0 * vector(2)
        result(2) = x(1) * vector(1) + 2.0d0 * x(2) * vector(2)
        res = 0
    end function jacobian_operator_hs15

    ! Jacobian transposed operator
    function jacobian_transposed_operator_hs15(number_variables, number_constraints, x, evaluate_at_x, &
                                               vector, result, user_data) result(res) &
        bind(C)
        integer(uno_int), value :: number_variables, number_constraints
        real(c_double), intent(in) :: x(*)
        logical(c_bool), value :: evaluate_at_x
        real(c_double), intent(in) :: vector(*)
        real(c_double), intent(out) :: result(*)
        type(c_ptr), value :: user_data
        integer(uno_int) :: res

        result(1) = x(2) * vector(1) + x(1) * vector(2)
        result(2) = 1.0d0 * vector(1) + 2.0d0 * x(2) * vector(2)
        res = 0
    end function jacobian_transposed_operator_hs15

    ! Lagrangian Hessian
    function lagrangian_hessian_hs15(number_variables, number_constraints, number_hessian_nonzeros, &
                                     x, objective_multiplier, multipliers, hessian_values, user_data) result(res) &
        bind(C)
        integer(uno_int), value :: number_variables, number_constraints, number_hessian_nonzeros
        real(c_double), intent(in) :: x(*), multipliers(*)
        real(c_double), intent(out) :: hessian_values(*)
        real(c_double), value :: objective_multiplier
        type(c_ptr), value :: user_data
        integer(uno_int) :: res

        hessian_values(1) = objective_multiplier * (1200.0d0 * x(1)**2 - 400.0d0 * x(2) + 2.0d0)
        hessian_values(2) = -400.0d0 * objective_multiplier * x(1) - multipliers(1)
        hessian_values(3) = 200.0d0 * objective_multiplier - 2.0d0 * multipliers(2)
        res = 0
    end function lagrangian_hessian_hs15

    ! Lagrangian Hessian operator
    function lagrangian_hessian_operator_hs15(number_variables, number_constraints, x, evaluate_at_x, &
                                              objective_multiplier, multipliers, vector, result, user_data) result(res) bind(C)
        integer(uno_int), value :: number_variables, number_constraints
        real(c_double), intent(in) :: x(*)
        logical(c_bool), value :: evaluate_at_x
        real(c_double), value :: objective_multiplier
        real(c_double), intent(in) :: multipliers(*), vector(*)
        real(c_double), intent(out) :: result(*)
        type(c_ptr), value :: user_data
        integer(uno_int) :: res
        real(c_double) :: hessian11, hessian12, hessian22

        hessian11 = objective_multiplier * (1200.0d0 * x(1)**2 - 400.0d0 * x(2) + 2.0d0)
        hessian12 = -400.0d0 * objective_multiplier * x(1) - multipliers(1)
        hessian22 = 200.0d0 * objective_multiplier - 2.0d0 * multipliers(2)
        result(1) = hessian11 * vector(1) + hessian12 * vector(2)
        result(2) = hessian12 * vector(1) + hessian22 * vector(2)
        res = 0
    end function lagrangian_hessian_operator_hs15

end program example_uno
