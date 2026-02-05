program example_uno
    include 'uno.f90'

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
    real(c_double), dimension(number_constraints) :: constraint_dual_solution
    real(c_double) :: solution_primal_feasibility, solution_stationarity, solution_complementarity
    integer(uno_int) :: optimization_status, iterate_status
    logical(c_bool) :: success, print_solution = .true.
    integer(uno_int) :: max_iterations = 1000
    real(c_double) :: primal_tolerance = 1.0d-6
    character(len=*), parameter :: hessian_model = "exact"//c_null_char
    integer(uno_int), parameter :: base_indexing = UNO_ONE_BASED_INDEXING
    integer(uno_int), parameter :: optimization_sense = UNO_MINIMIZE
    character(len=1), parameter :: hessian_triangular_part = UNO_LOWER_TRIANGLE
    real(c_double), parameter :: lagrangian_sign_convention = UNO_MULTIPLIER_NEGATIVE
    type(c_funptr) :: objective, gradient, constraints, jacobian, lagrangian_hessian
    type(c_funptr) :: jacobian_operator, jacobian_transposed_operator, lagrangian_hessian_operator

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
    objective = c_funloc(objective_hs15)
    gradient = c_funloc(gradient_hs15)
    constraints = c_funloc(constraints_hs15)
    jacobian = c_funloc(jacobian_hs15)
    jacobian_operator = c_funloc(jacobian_operator_hs15)
    jacobian_transposed_operator = c_funloc(jacobian_transposed_operator_hs15)
    lagrangian_hessian = c_funloc(lagrangian_hessian_hs15)
    lagrangian_hessian_operator = c_funloc(lagrangian_hessian_operator_hs15)

    !---------------------------------------------------
    ! Model creation
    !---------------------------------------------------
    model = uno_create_model(UNO_PROBLEM_NONLINEAR, number_variables, variables_lower_bounds, variables_upper_bounds, base_indexing)

    success = uno_set_objective(model, optimization_sense, objective, gradient)
    success = uno_set_constraints(model, number_constraints, constraints, constraints_lower_bounds, constraints_upper_bounds, &
                                  number_jacobian_nonzeros, jacobian_row_indices, jacobian_column_indices, jacobian)
    success = uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part, hessian_row_indices, &
                                         hessian_column_indices, lagrangian_hessian, lagrangian_sign_convention)
    success = uno_set_jacobian_operator(model, jacobian_operator)
    success = uno_set_jacobian_transposed_operator(model, jacobian_transposed_operator)
    success = uno_set_lagrangian_hessian_operator(model, lagrangian_hessian_operator, lagrangian_sign_convention)
    success = uno_set_initial_primal_iterate(model, x0)
    success = uno_set_initial_dual_iterate(model, y0)

    !---------------------------------------------------
    ! Solver creation
    !---------------------------------------------------
    solver = uno_create_solver()
    success = uno_set_solver_integer_option(solver, "max_iterations"//c_null_char, max_iterations)
    success = uno_set_solver_double_option(solver, "primal_tolerance"//c_null_char, primal_tolerance)
    success = uno_set_solver_bool_option(solver, "print_solution"//c_null_char, print_solution)
    success = uno_set_solver_string_option(solver, "hessian_model"//c_null_char, hessian_model)
    success = uno_set_solver_preset(solver, "filtersqp"//c_null_char)

    !---------------------------------------------------
    ! Solve
    !---------------------------------------------------
    call uno_optimize(solver, model)

    !---------------------------------------------------
    ! Get solution
    !---------------------------------------------------
    optimization_status = uno_get_optimization_status(solver)
    print *, 'optimization_status = ', optimization_status

    iterate_status = uno_get_solution_status(solver)
    print *, 'iterate_status = ', iterate_status

    solution_objective = uno_get_solution_objective(solver)
    print *, 'Solution objective = ', solution_objective

    call uno_get_primal_solution(solver, primal_solution)
    print *, 'Primal solution = ', primal_solution

    call uno_get_constraint_dual_solution(solver, constraint_dual_solution)
    print *, 'Constraint dual solution = ', constraint_dual_solution

    call uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)
    print *, 'Lower bound dual solution = ', lower_bound_dual_solution

    call uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)
    print *, 'Upper bound dual solution = ', upper_bound_dual_solution

    solution_primal_feasibility = uno_get_solution_primal_feasibility(solver)
    print *, 'Primal feasibility at solution = ', solution_primal_feasibility

    solution_stationarity = uno_get_solution_stationarity(solver)
    print *, 'Stationarity at solution = ', solution_stationarity

    solution_complementarity = uno_get_solution_complementarity(solver)
    print *, 'Complementarity at solution = ', solution_complementarity

    !---------------------------------------------------
    ! Cleanup
    !---------------------------------------------------
    call uno_destroy_solver(solver)
    call uno_destroy_model(model)

contains

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
        logical, value :: evaluate_at_x
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
        logical, value :: evaluate_at_x
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
        logical, value :: evaluate_at_x
        real(c_double), value :: objective_multiplier
        real(c_double), intent(in) :: multipliers(*), vector(*)
        real(c_double), intent(out) :: result(*)
        type(c_ptr), value :: user_data
        integer(uno_int) :: res
        real(c_double) :: hessian00, hessian10, hessian11

        hessian11 = objective_multiplier * (1200.0d0 * x(1)**2 - 400.0d0 * x(2) + 2.0d0)
        hessian12 = -400.0d0 * objective_multiplier * x(1) - multipliers(1)
        hessian22 = 200.0d0 * objective_multiplier - 2.0d0 * multipliers(2)
        result(1) = hessian11 * vector(1) + hessian12 * vector(2)
        result(2) = hessian12 * vector(1) + hessian22 * vector(2)
        res = 0
    end function lagrangian_hessian_operator_hs15

end program example_uno
