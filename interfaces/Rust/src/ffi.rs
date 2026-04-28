// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#![allow(non_camel_case_types)]

use std::os::raw::{c_char, c_int, c_void};

pub type uno_int = i32;

// ─────────────────────────────────────────────
// Optimization sense
// ─────────────────────────────────────────────
pub const UNO_MINIMIZE: uno_int =  1;
pub const UNO_MAXIMIZE: uno_int = -1;

// ─────────────────────────────────────────────
// Lagrange multiplier sign convention
// ─────────────────────────────────────────────
pub const UNO_MULTIPLIER_POSITIVE: uno_int =  1;
pub const UNO_MULTIPLIER_NEGATIVE: uno_int = -1;

// ─────────────────────────────────────────────
// Problem types
// ─────────────────────────────────────────────
pub const UNO_PROBLEM_LINEAR:    &[u8] = b"LP\0";
pub const UNO_PROBLEM_QUADRATIC: &[u8] = b"QP\0";
pub const UNO_PROBLEM_NONLINEAR: &[u8] = b"NLP\0";

// ─────────────────────────────────────────────
// Base indexing style
// ─────────────────────────────────────────────
pub const UNO_ZERO_BASED_INDEXING: uno_int = 0;
pub const UNO_ONE_BASED_INDEXING:  uno_int = 1;

// ─────────────────────────────────────────────
// Triangular part
// ─────────────────────────────────────────────
pub const UNO_LOWER_TRIANGLE: c_char = b'L' as c_char;
pub const UNO_UPPER_TRIANGLE: c_char = b'U' as c_char;

// ─────────────────────────────────────────────
// Option type
// ─────────────────────────────────────────────
pub const UNO_OPTION_TYPE_INTEGER:   uno_int =  0;
pub const UNO_OPTION_TYPE_DOUBLE:    uno_int =  1;
pub const UNO_OPTION_TYPE_BOOL:      uno_int =  2;
pub const UNO_OPTION_TYPE_STRING:    uno_int =  3;
pub const UNO_OPTION_TYPE_NOT_FOUND: uno_int = -1;

// ─────────────────────────────────────────────
// Optimization status
// ─────────────────────────────────────────────
pub const UNO_SUCCESS:           uno_int = 0;
pub const UNO_ITERATION_LIMIT:   uno_int = 1;
pub const UNO_TIME_LIMIT:        uno_int = 2;
pub const UNO_EVALUATION_ERROR:  uno_int = 3;
pub const UNO_ALGORITHMIC_ERROR: uno_int = 4;
pub const UNO_USER_TERMINATION:  uno_int = 5;

// ─────────────────────────────────────────────
// Solution status
// ─────────────────────────────────────────────
pub const UNO_NOT_OPTIMAL:                 uno_int = 0;
pub const UNO_FEASIBLE_KKT_POINT:          uno_int = 1;
pub const UNO_FEASIBLE_FJ_POINT:           uno_int = 2;
pub const UNO_INFEASIBLE_STATIONARY_POINT: uno_int = 3;
pub const UNO_FEASIBLE_SMALL_STEP:         uno_int = 4;
pub const UNO_INFEASIBLE_SMALL_STEP:       uno_int = 5;
pub const UNO_UNBOUNDED:                   uno_int = 6;

// ─────────────────────────────────────────────
// Version
// ─────────────────────────────────────────────
pub const UNO_VERSION_MAJOR: uno_int = 2;
pub const UNO_VERSION_MINOR: uno_int = 7;
pub const UNO_VERSION_PATCH: uno_int = 1;

// ─────────────────────────────────────────────
// Callback type aliases
// ─────────────────────────────────────────────

pub type uno_objective_callback = unsafe extern "C" fn(
    number_variables: uno_int,
    x: *const f64,
    objective_value: *mut f64,
    user_data: *mut c_void,
) -> uno_int;

pub type uno_constraints_callback = unsafe extern "C" fn(
    number_variables: uno_int,
    number_constraints: uno_int,
    x: *const f64,
    constraint_values: *mut f64,
    user_data: *mut c_void,
) -> uno_int;

pub type uno_objective_gradient_callback = unsafe extern "C" fn(
    number_variables: uno_int,
    x: *const f64,
    gradient: *mut f64,
    user_data: *mut c_void,
) -> uno_int;

pub type uno_constraints_jacobian_callback = unsafe extern "C" fn(
    number_variables: uno_int,
    number_jacobian_nonzeros: uno_int,
    x: *const f64,
    jacobian_values: *mut f64,
    user_data: *mut c_void,
) -> uno_int;

pub type uno_lagrangian_hessian_callback = unsafe extern "C" fn(
    number_variables: uno_int,
    number_constraints: uno_int,
    number_hessian_nonzeros: uno_int,
    x: *const f64,
    objective_multiplier: f64,
    multipliers: *const f64,
    hessian_values: *mut f64,
    user_data: *mut c_void,
) -> uno_int;

pub type uno_constraints_jacobian_operator_callback = unsafe extern "C" fn(
    number_variables: uno_int,
    number_constraints: uno_int,
    x: *const f64,
    evaluate_at_x: bool,
    vector: *const f64,
    result: *mut f64,
    user_data: *mut c_void,
) -> uno_int;

pub type uno_constraints_jacobian_transposed_operator_callback = unsafe extern "C" fn(
    number_variables: uno_int,
    number_constraints: uno_int,
    x: *const f64,
    evaluate_at_x: bool,
    vector: *const f64,
    result: *mut f64,
    user_data: *mut c_void,
) -> uno_int;

pub type uno_lagrangian_hessian_operator_callback = unsafe extern "C" fn(
    number_variables: uno_int,
    number_constraints: uno_int,
    x: *const f64,
    evaluate_at_x: bool,
    objective_multiplier: f64,
    multipliers: *const f64,
    vector: *const f64,
    result: *mut f64,
    user_data: *mut c_void,
) -> uno_int;

pub type uno_notify_acceptable_iterate_callback = unsafe extern "C" fn(
    number_variables: uno_int,
    number_constraints: uno_int,
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

pub type uno_termination_callback = unsafe extern "C" fn(
    number_variables: uno_int,
    number_constraints: uno_int,
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

pub type uno_logger_stream_callback = unsafe extern "C" fn(
    buffer: *const c_char,
    length: uno_int,
    user_data: *mut c_void,
) -> uno_int;

// ─────────────────────────────────────────────
// extern "C" declarations
// ─────────────────────────────────────────────

#[link(name = "uno")]
extern "C" {
    // uno_get_version
    pub fn uno_get_version(major: *mut uno_int, minor: *mut uno_int, patch: *mut uno_int);

    // uno_create_model
    pub fn uno_create_model(
        problem_type: *const c_char,
        number_variables: uno_int,
        variables_lower_bounds: *const f64,
        variables_upper_bounds: *const f64,
        base_indexing: uno_int,
    ) -> *mut c_void;

    // uno_set_objective
    pub fn uno_set_objective(
        model: *mut c_void,
        optimization_sense: uno_int,
        objective_function: uno_objective_callback,
        objective_gradient: uno_objective_gradient_callback,
    ) -> bool;

    // uno_set_constraints
    pub fn uno_set_constraints(
        model: *mut c_void,
        number_constraints: uno_int,
        constraint_functions: uno_constraints_callback,
        constraints_lower_bounds: *const f64,
        constraints_upper_bounds: *const f64,
        number_jacobian_nonzeros: uno_int,
        jacobian_row_indices: *const uno_int,
        jacobian_column_indices: *const uno_int,
        jacobian: uno_constraints_jacobian_callback,
    ) -> bool;

    // uno_set_jacobian_operator
    pub fn uno_set_jacobian_operator(
        model: *mut c_void,
        jacobian_operator: uno_constraints_jacobian_operator_callback,
    ) -> bool;

    // uno_set_jacobian_transposed_operator
    pub fn uno_set_jacobian_transposed_operator(
        model: *mut c_void,
        jacobian_transposed_operator: uno_constraints_jacobian_transposed_operator_callback,
    ) -> bool;

    // uno_set_lagrangian_hessian
    pub fn uno_set_lagrangian_hessian(
        model: *mut c_void,
        number_hessian_nonzeros: uno_int,
        hessian_triangular_part: c_char,
        hessian_row_indices: *const uno_int,
        hessian_column_indices: *const uno_int,
        lagrangian_hessian: uno_lagrangian_hessian_callback,
    ) -> bool;

    // uno_set_lagrangian_hessian_operator
    pub fn uno_set_lagrangian_hessian_operator(
        model: *mut c_void,
        lagrangian_hessian_operator: uno_lagrangian_hessian_operator_callback,
    ) -> bool;

    // uno_set_lagrangian_sign_convention
    pub fn uno_set_lagrangian_sign_convention(
        model: *mut c_void,
        lagrangian_sign_convention: uno_int,
    ) -> bool;

    // uno_set_user_data
    pub fn uno_set_user_data(model: *mut c_void, user_data: *mut c_void) -> bool;

    // uno_set_initial_primal_iterate_component
    pub fn uno_set_initial_primal_iterate_component(
        model: *mut c_void,
        index: uno_int,
        initial_primal_component: f64,
    ) -> bool;

    // uno_set_initial_dual_iterate_component
    pub fn uno_set_initial_dual_iterate_component(
        model: *mut c_void,
        index: uno_int,
        initial_dual_component: f64,
    ) -> bool;

    // uno_set_initial_primal_iterate
    pub fn uno_set_initial_primal_iterate(
        model: *mut c_void,
        initial_primal_iterate: *const f64,
    ) -> bool;

    // uno_set_initial_dual_iterate
    pub fn uno_set_initial_dual_iterate(
        model: *mut c_void,
        initial_dual_iterate: *const f64,
    ) -> bool;

    // uno_create_solver
    pub fn uno_create_solver() -> *mut c_void;

    // uno_set_solver_integer_option
    pub fn uno_set_solver_integer_option(
        solver: *mut c_void,
        option_name: *const c_char,
        option_value: uno_int,
    ) -> bool;

    // uno_set_solver_double_option
    pub fn uno_set_solver_double_option(
        solver: *mut c_void,
        option_name: *const c_char,
        option_value: f64,
    ) -> bool;

    // uno_set_solver_bool_option
    pub fn uno_set_solver_bool_option(
        solver: *mut c_void,
        option_name: *const c_char,
        option_value: bool,
    ) -> bool;

    // uno_set_solver_string_option
    pub fn uno_set_solver_string_option(
        solver: *mut c_void,
        option_name: *const c_char,
        option_value: *const c_char,
    ) -> bool;

    // uno_get_solver_option_type
    pub fn uno_get_solver_option_type(solver: *mut c_void, option_name: *const c_char) -> uno_int;

    // uno_load_solver_option_file
    pub fn uno_load_solver_option_file(solver: *mut c_void, file_name: *const c_char) -> bool;

    // uno_set_solver_preset
    pub fn uno_set_solver_preset(solver: *mut c_void, preset_name: *const c_char) -> bool;

    // uno_set_solver_callbacks
    pub fn uno_set_solver_callbacks(
        solver: *mut c_void,
        notify_acceptable_iterate_callback: Option<uno_notify_acceptable_iterate_callback>,
        termination_callback: Option<uno_termination_callback>,
        user_data: *mut c_void,
    ) -> bool;

    // uno_set_logger_stream_callback
    pub fn uno_set_logger_stream_callback(
        logger_stream_callback: uno_logger_stream_callback,
        user_data: *mut c_void,
    ) -> bool;

    // uno_reset_logger_stream
    pub fn uno_reset_logger_stream() -> bool;

    // uno_optimize
    pub fn uno_optimize(solver: *mut c_void, model: *mut c_void);

    // uno_get_solver_integer_option
    pub fn uno_get_solver_integer_option(solver: *mut c_void, option_name: *const c_char) -> c_int;

    // uno_get_solver_double_option
    pub fn uno_get_solver_double_option(solver: *mut c_void, option_name: *const c_char) -> f64;

    // uno_get_solver_bool_option
    pub fn uno_get_solver_bool_option(solver: *mut c_void, option_name: *const c_char) -> bool;

    // uno_get_solver_string_option
    pub fn uno_get_solver_string_option(
        solver: *mut c_void,
        option_name: *const c_char,
    ) -> *const c_char;

    // uno_get_optimization_status
    pub fn uno_get_optimization_status(solver: *mut c_void) -> uno_int;

    // uno_get_solution_status
    pub fn uno_get_solution_status(solver: *mut c_void) -> uno_int;

    // uno_get_solution_objective
    pub fn uno_get_solution_objective(solver: *mut c_void) -> f64;

    // uno_get_primal_solution_component
    pub fn uno_get_primal_solution_component(solver: *mut c_void, index: uno_int) -> f64;

    // uno_get_constraint_dual_solution_component
    pub fn uno_get_constraint_dual_solution_component(solver: *mut c_void, index: uno_int) -> f64;

    // uno_get_lower_bound_dual_solution_component
    pub fn uno_get_lower_bound_dual_solution_component(solver: *mut c_void, index: uno_int) -> f64;

    // uno_get_upper_bound_dual_solution_component
    pub fn uno_get_upper_bound_dual_solution_component(solver: *mut c_void, index: uno_int) -> f64;

    // uno_get_primal_solution
    pub fn uno_get_primal_solution(solver: *mut c_void, primal_solution: *mut f64);

    // uno_get_constraint_dual_solution
    pub fn uno_get_constraint_dual_solution(
        solver: *mut c_void,
        constraint_dual_solution: *mut f64,
    );

    // uno_get_lower_bound_dual_solution
    pub fn uno_get_lower_bound_dual_solution(
        solver: *mut c_void,
        lower_bound_dual_solution: *mut f64,
    );

    // uno_get_upper_bound_dual_solution
    pub fn uno_get_upper_bound_dual_solution(
        solver: *mut c_void,
        upper_bound_dual_solution: *mut f64,
    );

    // uno_get_solution_primal_feasibility
    pub fn uno_get_solution_primal_feasibility(solver: *mut c_void) -> f64;

    // uno_get_solution_stationarity
    pub fn uno_get_solution_stationarity(solver: *mut c_void) -> f64;

    // uno_get_solution_complementarity
    pub fn uno_get_solution_complementarity(solver: *mut c_void) -> f64;

    // uno_get_number_iterations
    pub fn uno_get_number_iterations(solver: *mut c_void) -> uno_int;

    // uno_get_cpu_time
    pub fn uno_get_cpu_time(solver: *mut c_void) -> f64;

    // uno_get_number_objective_evaluations
    pub fn uno_get_number_objective_evaluations(solver: *mut c_void) -> uno_int;

    // uno_get_number_constraint_evaluations
    pub fn uno_get_number_constraint_evaluations(solver: *mut c_void) -> uno_int;

    // uno_get_number_objective_gradient_evaluations
    pub fn uno_get_number_objective_gradient_evaluations(solver: *mut c_void) -> uno_int;

    // uno_get_number_jacobian_evaluations
    pub fn uno_get_number_jacobian_evaluations(solver: *mut c_void) -> uno_int;

    // uno_get_number_hessian_evaluations
    pub fn uno_get_number_hessian_evaluations(solver: *mut c_void) -> uno_int;

    // uno_get_number_subproblem_solved_evaluations
    pub fn uno_get_number_subproblem_solved_evaluations(solver: *mut c_void) -> uno_int;

    // uno_destroy_model
    pub fn uno_destroy_model(model: *mut c_void);

    // uno_destroy_solver
    pub fn uno_destroy_solver(solver: *mut c_void);
}
