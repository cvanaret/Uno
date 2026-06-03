// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#![allow(non_camel_case_types)]

use std::os::raw::{c_char, c_int, c_void};

//--- UNO_INT_TYPE ---

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
