// Raw FFI bindings for the Uno C API
// Mirrors Uno_C_API.h

use std::os::raw::{c_char, c_int, c_void};

pub type UnoInt = i32;

// Optimization sense
pub const UNO_MINIMIZE: UnoInt = 1;
pub const UNO_MAXIMIZE: UnoInt = -1;

// Lagrange multiplier sign convention
pub const UNO_MULTIPLIER_POSITIVE: UnoInt = 1;
pub const UNO_MULTIPLIER_NEGATIVE: UnoInt = -1;

// Problem types
pub const UNO_PROBLEM_NONLINEAR: &[u8] = b"NLP\0";
pub const UNO_PROBLEM_QUADRATIC: &[u8] = b"QP\0";
pub const UNO_PROBLEM_LINEAR: &[u8] = b"LP\0";

// Base indexing
pub const UNO_ZERO_BASED_INDEXING: UnoInt = 0;
pub const UNO_ONE_BASED_INDEXING: UnoInt = 1;

// Triangular part
pub const UNO_LOWER_TRIANGLE: c_char = b'L' as c_char;
pub const UNO_UPPER_TRIANGLE: c_char = b'U' as c_char;

// Optimization status
pub const UNO_SUCCESS: UnoInt = 0;
pub const UNO_ITERATION_LIMIT: UnoInt = 1;
pub const UNO_TIME_LIMIT: UnoInt = 2;
pub const UNO_EVALUATION_ERROR: UnoInt = 3;
pub const UNO_ALGORITHMIC_ERROR: UnoInt = 4;
pub const UNO_USER_TERMINATION: UnoInt = 5;

// Iterate (solution) status
pub const UNO_NOT_OPTIMAL: UnoInt = 0;
pub const UNO_FEASIBLE_KKT_POINT: UnoInt = 1;
pub const UNO_FEASIBLE_FJ_POINT: UnoInt = 2;
pub const UNO_INFEASIBLE_STATIONARY_POINT: UnoInt = 3;
pub const UNO_FEASIBLE_SMALL_STEP: UnoInt = 4;
pub const UNO_INFEASIBLE_SMALL_STEP: UnoInt = 5;
pub const UNO_UNBOUNDED: UnoInt = 6;

// ──────────────────────────────────────────────────────────────────────────────
// Callback type aliases
// ──────────────────────────────────────────────────────────────────────────────

pub type UnoObjectiveCallback = unsafe extern "C" fn(
    number_variables: UnoInt,
    x: *const f64,
    objective_value: *mut f64,
    user_data: *mut c_void,
) -> UnoInt;

pub type UnoConstraintsCallback = unsafe extern "C" fn(
    number_variables: UnoInt,
    number_constraints: UnoInt,
    x: *const f64,
    constraint_values: *mut f64,
    user_data: *mut c_void,
) -> UnoInt;

pub type UnoObjectiveGradientCallback = unsafe extern "C" fn(
    number_variables: UnoInt,
    x: *const f64,
    gradient: *mut f64,
    user_data: *mut c_void,
) -> UnoInt;

pub type UnoConstraintsJacobianCallback = unsafe extern "C" fn(
    number_variables: UnoInt,
    number_jacobian_nonzeros: UnoInt,
    x: *const f64,
    jacobian_values: *mut f64,
    user_data: *mut c_void,
) -> UnoInt;

pub type UnoLagrangianHessianCallback = unsafe extern "C" fn(
    number_variables: UnoInt,
    number_constraints: UnoInt,
    number_hessian_nonzeros: UnoInt,
    x: *const f64,
    objective_multiplier: f64,
    multipliers: *const f64,
    hessian_values: *mut f64,
    user_data: *mut c_void,
) -> UnoInt;

pub type UnoConstraintsJacobianOperatorCallback = unsafe extern "C" fn(
    number_variables: UnoInt,
    number_constraints: UnoInt,
    x: *const f64,
    evaluate_at_x: bool,
    vector: *const f64,
    result: *mut f64,
    user_data: *mut c_void,
) -> UnoInt;

pub type UnoConstraintsJacobianTransposedOperatorCallback = unsafe extern "C" fn(
    number_variables: UnoInt,
    number_constraints: UnoInt,
    x: *const f64,
    evaluate_at_x: bool,
    vector: *const f64,
    result: *mut f64,
    user_data: *mut c_void,
) -> UnoInt;

pub type UnoLagrangianHessianOperatorCallback = unsafe extern "C" fn(
    number_variables: UnoInt,
    number_constraints: UnoInt,
    x: *const f64,
    evaluate_at_x: bool,
    objective_multiplier: f64,
    multipliers: *const f64,
    vector: *const f64,
    result: *mut f64,
    user_data: *mut c_void,
) -> UnoInt;

pub type UnoNotifyAcceptableIterateCallback = unsafe extern "C" fn(
    number_variables: UnoInt,
    number_constraints: UnoInt,
    primals: *const f64,
    lower_bound_multipliers: *const f64,
    upper_bound_multipliers: *const f64,
    constraint_multipliers: *const f64,
    objective_multiplier: f64,
    primal_feasibility_residual: f64,
    stationarity_residual: f64,
    complementarity_residual: f64,
    user_data: *mut c_void,
);

pub type UnoTerminationCallback = unsafe extern "C" fn(
    number_variables: UnoInt,
    number_constraints: UnoInt,
    primals: *const f64,
    lower_bound_multipliers: *const f64,
    upper_bound_multipliers: *const f64,
    constraint_multipliers: *const f64,
    objective_multiplier: f64,
    primal_feasibility_residual: f64,
    stationarity_residual: f64,
    complementarity_residual: f64,
    user_data: *mut c_void,
) -> bool;

pub type UnoLoggerStreamCallback = unsafe extern "C" fn(
    buffer: *const c_char,
    length: UnoInt,
    user_data: *mut c_void,
) -> UnoInt;

// ──────────────────────────────────────────────────────────────────────────────
// Extern "C" declarations
// ──────────────────────────────────────────────────────────────────────────────

#[link(name = "uno")]
extern "C" {
    // Version
    pub fn uno_get_version(major: *mut UnoInt, minor: *mut UnoInt, patch: *mut UnoInt);

    // Model lifecycle
    pub fn uno_create_model(
        problem_type: *const c_char,
        number_variables: UnoInt,
        variables_lower_bounds: *const f64,
        variables_upper_bounds: *const f64,
        base_indexing: UnoInt,
    ) -> *mut c_void;

    pub fn uno_destroy_model(model: *mut c_void);

    // Model setup
    pub fn uno_set_objective(
        model: *mut c_void,
        optimization_sense: UnoInt,
        objective_function: UnoObjectiveCallback,
        objective_gradient: UnoObjectiveGradientCallback,
    ) -> bool;

    pub fn uno_set_constraints(
        model: *mut c_void,
        number_constraints: UnoInt,
        constraint_functions: UnoConstraintsCallback,
        constraints_lower_bounds: *const f64,
        constraints_upper_bounds: *const f64,
        number_jacobian_nonzeros: UnoInt,
        jacobian_row_indices: *const UnoInt,
        jacobian_column_indices: *const UnoInt,
        jacobian: UnoConstraintsJacobianCallback,
    ) -> bool;

    pub fn uno_set_jacobian_operator(
        model: *mut c_void,
        jacobian_operator: UnoConstraintsJacobianOperatorCallback,
    ) -> bool;

    pub fn uno_set_jacobian_transposed_operator(
        model: *mut c_void,
        jacobian_transposed_operator: UnoConstraintsJacobianTransposedOperatorCallback,
    ) -> bool;

    pub fn uno_set_lagrangian_hessian(
        model: *mut c_void,
        number_hessian_nonzeros: UnoInt,
        hessian_triangular_part: c_char,
        hessian_row_indices: *const UnoInt,
        hessian_column_indices: *const UnoInt,
        lagrangian_hessian: UnoLagrangianHessianCallback,
    ) -> bool;

    pub fn uno_set_lagrangian_hessian_operator(
        model: *mut c_void,
        lagrangian_hessian_operator: UnoLagrangianHessianOperatorCallback,
    ) -> bool;

    pub fn uno_set_lagrangian_sign_convention(
        model: *mut c_void,
        lagrangian_sign_convention: UnoInt,
    ) -> bool;

    pub fn uno_set_user_data(model: *mut c_void, user_data: *mut c_void) -> bool;

    pub fn uno_set_initial_primal_iterate(
        model: *mut c_void,
        initial_primal_iterate: *const f64,
    ) -> bool;

    pub fn uno_set_initial_dual_iterate(
        model: *mut c_void,
        initial_dual_iterate: *const f64,
    ) -> bool;

    pub fn uno_set_initial_primal_iterate_component(
        model: *mut c_void,
        index: UnoInt,
        initial_primal_component: f64,
    ) -> bool;

    pub fn uno_set_initial_dual_iterate_component(
        model: *mut c_void,
        index: UnoInt,
        initial_dual_component: f64,
    ) -> bool;

    // Solver lifecycle
    pub fn uno_create_solver() -> *mut c_void;
    pub fn uno_destroy_solver(solver: *mut c_void);

    // Solver options
    pub fn uno_set_solver_integer_option(
        solver: *mut c_void,
        option_name: *const c_char,
        option_value: UnoInt,
    ) -> bool;
    pub fn uno_set_solver_double_option(
        solver: *mut c_void,
        option_name: *const c_char,
        option_value: f64,
    ) -> bool;
    pub fn uno_set_solver_bool_option(
        solver: *mut c_void,
        option_name: *const c_char,
        option_value: bool,
    ) -> bool;
    pub fn uno_set_solver_string_option(
        solver: *mut c_void,
        option_name: *const c_char,
        option_value: *const c_char,
    ) -> bool;

    pub fn uno_get_solver_option_type(solver: *mut c_void, option_name: *const c_char) -> UnoInt;
    pub fn uno_get_solver_integer_option(solver: *mut c_void, option_name: *const c_char) -> c_int;
    pub fn uno_get_solver_double_option(solver: *mut c_void, option_name: *const c_char) -> f64;
    pub fn uno_get_solver_bool_option(solver: *mut c_void, option_name: *const c_char) -> bool;
    pub fn uno_get_solver_string_option(
        solver: *mut c_void,
        option_name: *const c_char,
    ) -> *const c_char;

    pub fn uno_load_solver_option_file(solver: *mut c_void, file_name: *const c_char) -> bool;
    pub fn uno_set_solver_preset(solver: *mut c_void, preset_name: *const c_char) -> bool;
    pub fn uno_set_solver_callbacks(
        solver: *mut c_void,
        notify_acceptable_iterate_callback: Option<UnoNotifyAcceptableIterateCallback>,
        termination_callback: Option<UnoTerminationCallback>,
        user_data: *mut c_void,
    ) -> bool;

    pub fn uno_set_logger_stream_callback(
        logger_stream_callback: UnoLoggerStreamCallback,
        user_data: *mut c_void,
    ) -> bool;
    pub fn uno_reset_logger_stream() -> bool;

    // Solve
    pub fn uno_optimize(solver: *mut c_void, model: *mut c_void);

    // Results
    pub fn uno_get_optimization_status(solver: *mut c_void) -> UnoInt;
    pub fn uno_get_solution_status(solver: *mut c_void) -> UnoInt;
    pub fn uno_get_solution_objective(solver: *mut c_void) -> f64;

    pub fn uno_get_primal_solution_component(solver: *mut c_void, index: UnoInt) -> f64;
    pub fn uno_get_constraint_dual_solution_component(solver: *mut c_void, index: UnoInt) -> f64;
    pub fn uno_get_lower_bound_dual_solution_component(solver: *mut c_void, index: UnoInt) -> f64;
    pub fn uno_get_upper_bound_dual_solution_component(solver: *mut c_void, index: UnoInt) -> f64;

    pub fn uno_get_primal_solution(solver: *mut c_void, primal_solution: *mut f64);
    pub fn uno_get_constraint_dual_solution(
        solver: *mut c_void,
        constraint_dual_solution: *mut f64,
    );
    pub fn uno_get_lower_bound_dual_solution(
        solver: *mut c_void,
        lower_bound_dual_solution: *mut f64,
    );
    pub fn uno_get_upper_bound_dual_solution(
        solver: *mut c_void,
        upper_bound_dual_solution: *mut f64,
    );

    pub fn uno_get_solution_primal_feasibility(solver: *mut c_void) -> f64;
    pub fn uno_get_solution_stationarity(solver: *mut c_void) -> f64;
    pub fn uno_get_solution_complementarity(solver: *mut c_void) -> f64;

    pub fn uno_get_number_iterations(solver: *mut c_void) -> UnoInt;
    pub fn uno_get_cpu_time(solver: *mut c_void) -> f64;
    pub fn uno_get_number_objective_evaluations(solver: *mut c_void) -> UnoInt;
    pub fn uno_get_number_constraint_evaluations(solver: *mut c_void) -> UnoInt;
    pub fn uno_get_number_objective_gradient_evaluations(solver: *mut c_void) -> UnoInt;
    pub fn uno_get_number_jacobian_evaluations(solver: *mut c_void) -> UnoInt;
    pub fn uno_get_number_hessian_evaluations(solver: *mut c_void) -> UnoInt;
    pub fn uno_get_number_subproblem_solved_evaluations(solver: *mut c_void) -> UnoInt;
}
