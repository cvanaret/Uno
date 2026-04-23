// Rust wrapper over the Uno C API

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_void};

use crate::ffi;

// re-export status constants

pub use ffi::{
    UnoInt,
    UNO_MINIMIZE, UNO_MAXIMIZE,
    UNO_MULTIPLIER_POSITIVE, UNO_MULTIPLIER_NEGATIVE,
    UNO_ZERO_BASED_INDEXING, UNO_ONE_BASED_INDEXING,
    UNO_LOWER_TRIANGLE, UNO_UPPER_TRIANGLE,
    UNO_SUCCESS, UNO_ITERATION_LIMIT, UNO_TIME_LIMIT,
    UNO_EVALUATION_ERROR, UNO_ALGORITHMIC_ERROR, UNO_USER_TERMINATION,
    UNO_NOT_OPTIMAL, UNO_FEASIBLE_KKT_POINT, UNO_FEASIBLE_FJ_POINT,
    UNO_INFEASIBLE_STATIONARY_POINT, UNO_FEASIBLE_SMALL_STEP,
    UNO_INFEASIBLE_SMALL_STEP, UNO_UNBOUNDED,
};

// version

pub fn get_version() -> (i32, i32, i32) {
    let mut major: UnoInt = 0;
    let mut minor: UnoInt = 0;
    let mut patch: UnoInt = 0;
    unsafe { ffi::uno_get_version(&mut major, &mut minor, &mut patch) }
    (major, minor, patch)
}

// callback type re-exports

pub use ffi::{
    UnoObjectiveCallback,
    UnoObjectiveGradientCallback,
    UnoConstraintsCallback,
    UnoConstraintsJacobianCallback,
    UnoLagrangianHessianCallback,
    UnoConstraintsJacobianOperatorCallback,
    UnoConstraintsJacobianTransposedOperatorCallback,
    UnoLagrangianHessianOperatorCallback,
    UnoNotifyAcceptableIterateCallback,
    UnoTerminationCallback,
    UnoLoggerStreamCallback,
};

// Model

/// An optimization model that can be passed to `Solver::optimize`.
///
/// Owns the raw pointer returned by `uno_create_model` and destroys it on drop.
pub struct Model {
    ptr: *mut c_void,
    pub number_variables: UnoInt,
}

impl Drop for Model {
    fn drop(&mut self) {
        unsafe { ffi::uno_destroy_model(self.ptr) }
    }
}

impl Model {
    /// Creates a new model.
    ///
    /// # Arguments
    /// * `problem_type` – one of `"NLP"`, `"QP"`, or `"LP"`.
    /// * `lower_bounds` / `upper_bounds` – per-variable bounds (`f64::NEG_INFINITY` / `f64::INFINITY` for unbounded).
    /// * `base_indexing` – `UNO_ZERO_BASED_INDEXING` or `UNO_ONE_BASED_INDEXING`.
    pub fn new(
        problem_type: &str,
        lower_bounds: &[f64],
        upper_bounds: &[f64],
        base_indexing: UnoInt,
    ) -> Self {
        assert_eq!(
            lower_bounds.len(),
            upper_bounds.len(),
            "lower and upper bound arrays must have the same length"
        );
        let n = lower_bounds.len() as UnoInt;
        let pt = CString::new(problem_type).unwrap();
        let ptr = unsafe {
            ffi::uno_create_model(
                pt.as_ptr(),
                n,
                lower_bounds.as_ptr(),
                upper_bounds.as_ptr(),
                base_indexing,
            )
        };
        assert!(!ptr.is_null(), "uno_create_model returned null");
        Self { ptr, number_variables: n }
    }

    // objective

    pub fn set_objective(
        &mut self,
        sense: UnoInt,
        objective_fn: UnoObjectiveCallback,
        gradient_fn: UnoObjectiveGradientCallback,
    ) -> bool {
        unsafe { ffi::uno_set_objective(self.ptr, sense, objective_fn, gradient_fn) }
    }

    // constraints

    pub fn set_constraints(
        &mut self,
        constraint_fn: UnoConstraintsCallback,
        lower_bounds: &[f64],
        upper_bounds: &[f64],
        jacobian_row_indices: &[UnoInt],
        jacobian_col_indices: &[UnoInt],
        jacobian_fn: UnoConstraintsJacobianCallback,
    ) -> bool {
        let nc = lower_bounds.len() as UnoInt;
        let nnz = jacobian_row_indices.len() as UnoInt;
        assert_eq!(
            jacobian_row_indices.len(),
            jacobian_col_indices.len(),
            "jacobian row/col index arrays must have the same length"
        );
        unsafe {
            ffi::uno_set_constraints(
                self.ptr,
                nc,
                constraint_fn,
                lower_bounds.as_ptr(),
                upper_bounds.as_ptr(),
                nnz,
                jacobian_row_indices.as_ptr(),
                jacobian_col_indices.as_ptr(),
                jacobian_fn,
            )
        }
    }

    pub fn set_jacobian_operator(&mut self, op: UnoConstraintsJacobianOperatorCallback) -> bool {
        unsafe { ffi::uno_set_jacobian_operator(self.ptr, op) }
    }

    pub fn set_jacobian_transposed_operator(
        &mut self,
        op: UnoConstraintsJacobianTransposedOperatorCallback,
    ) -> bool {
        unsafe { ffi::uno_set_jacobian_transposed_operator(self.ptr, op) }
    }

    // Lagrangian Hessian

    pub fn set_lagrangian_hessian(
        &mut self,
        triangular_part: c_char,
        row_indices: &[UnoInt],
        col_indices: &[UnoInt],
        hessian_fn: UnoLagrangianHessianCallback,
    ) -> bool {
        let nnz = row_indices.len() as UnoInt;
        assert_eq!(row_indices.len(), col_indices.len());
        unsafe {
            ffi::uno_set_lagrangian_hessian(
                self.ptr,
                nnz,
                triangular_part,
                row_indices.as_ptr(),
                col_indices.as_ptr(),
                hessian_fn,
            )
        }
    }

    pub fn set_lagrangian_hessian_operator(
        &mut self,
        op: UnoLagrangianHessianOperatorCallback,
    ) -> bool {
        unsafe { ffi::uno_set_lagrangian_hessian_operator(self.ptr, op) }
    }

    pub fn set_lagrangian_sign_convention(&mut self, convention: UnoInt) -> bool {
        unsafe { ffi::uno_set_lagrangian_sign_convention(self.ptr, convention) }
    }

    // initial iterates

    pub fn set_initial_primal_iterate(&mut self, x0: &[f64]) -> bool {
        unsafe { ffi::uno_set_initial_primal_iterate(self.ptr, x0.as_ptr()) }
    }

    pub fn set_initial_dual_iterate(&mut self, y0: &[f64]) -> bool {
        unsafe { ffi::uno_set_initial_dual_iterate(self.ptr, y0.as_ptr()) }
    }

    pub fn set_initial_primal_component(&mut self, index: UnoInt, value: f64) -> bool {
        unsafe { ffi::uno_set_initial_primal_iterate_component(self.ptr, index, value) }
    }

    pub fn set_initial_dual_component(&mut self, index: UnoInt, value: f64) -> bool {
        unsafe { ffi::uno_set_initial_dual_iterate_component(self.ptr, index, value) }
    }

    // user data

    pub fn set_user_data(&mut self, data: *mut c_void) -> bool {
        unsafe { ffi::uno_set_user_data(self.ptr, data) }
    }
}

// Solution

/// The result returned by `Solver::optimize`.
pub struct Solution {
    pub optimization_status: UnoInt,
    pub solution_status: UnoInt,
    pub objective: f64,
    pub primals: Vec<f64>,
    pub constraint_duals: Vec<f64>,
    pub lower_bound_duals: Vec<f64>,
    pub upper_bound_duals: Vec<f64>,
    pub primal_feasibility: f64,
    pub stationarity: f64,
    pub complementarity: f64,
    pub iterations: UnoInt,
    pub cpu_time: f64,
}

impl Solution {
    pub fn optimization_status_str(&self) -> &'static str {
        match self.optimization_status {
            UNO_SUCCESS => "SUCCESS",
            UNO_ITERATION_LIMIT => "ITERATION_LIMIT",
            UNO_TIME_LIMIT => "TIME_LIMIT",
            UNO_EVALUATION_ERROR => "EVALUATION_ERROR",
            UNO_ALGORITHMIC_ERROR => "ALGORITHMIC_ERROR",
            UNO_USER_TERMINATION => "USER_TERMINATION",
            _ => "UNKNOWN",
        }
    }

    pub fn solution_status_str(&self) -> &'static str {
        match self.solution_status {
            UNO_NOT_OPTIMAL => "NOT_OPTIMAL",
            UNO_FEASIBLE_KKT_POINT => "FEASIBLE_KKT_POINT",
            UNO_FEASIBLE_FJ_POINT => "FEASIBLE_FJ_POINT",
            UNO_INFEASIBLE_STATIONARY_POINT => "INFEASIBLE_STATIONARY_POINT",
            UNO_FEASIBLE_SMALL_STEP => "FEASIBLE_SMALL_STEP",
            UNO_INFEASIBLE_SMALL_STEP => "INFEASIBLE_SMALL_STEP",
            UNO_UNBOUNDED => "UNBOUNDED",
            _ => "UNKNOWN",
        }
    }
}

// Solver

/// Wraps the Uno solver pointer.
pub struct Solver {
    ptr: *mut c_void,
}

impl Drop for Solver {
    fn drop(&mut self) {
        unsafe { ffi::uno_destroy_solver(self.ptr) }
    }
}

impl Solver {
    /// Creates a new solver.
    pub fn new() -> Self {
        let ptr = unsafe { ffi::uno_create_solver() };
        assert!(!ptr.is_null(), "uno_create_solver returned null");
        Self { ptr }
    }

    // option setters

    pub fn set_integer_option(&self, name: &str, value: UnoInt) -> bool {
        let name = CString::new(name).unwrap();
        unsafe { ffi::uno_set_solver_integer_option(self.ptr, name.as_ptr(), value) }
    }

    pub fn set_double_option(&self, name: &str, value: f64) -> bool {
        let name = CString::new(name).unwrap();
        unsafe { ffi::uno_set_solver_double_option(self.ptr, name.as_ptr(), value) }
    }

    pub fn set_bool_option(&self, name: &str, value: bool) -> bool {
        let name = CString::new(name).unwrap();
        unsafe { ffi::uno_set_solver_bool_option(self.ptr, name.as_ptr(), value) }
    }

    pub fn set_string_option(&self, name: &str, value: &str) -> bool {
        let name = CString::new(name).unwrap();
        let value = CString::new(value).unwrap();
        unsafe { ffi::uno_set_solver_string_option(self.ptr, name.as_ptr(), value.as_ptr()) }
    }

    pub fn set_preset(&self, preset: &str) -> bool {
        let preset = CString::new(preset).unwrap();
        unsafe { ffi::uno_set_solver_preset(self.ptr, preset.as_ptr()) }
    }

    pub fn load_option_file(&self, path: &str) -> bool {
        let path = CString::new(path).unwrap();
        unsafe { ffi::uno_load_solver_option_file(self.ptr, path.as_ptr()) }
    }

    // option getters

    pub fn get_string_option(&self, name: &str) -> String {
        let name = CString::new(name).unwrap();
        let ptr = unsafe { ffi::uno_get_solver_string_option(self.ptr, name.as_ptr()) };
        if ptr.is_null() {
            String::new()
        } else {
            unsafe { CStr::from_ptr(ptr).to_string_lossy().into_owned() }
        }
    }

    pub fn get_bool_option(&self, name: &str) -> bool {
        let name = CString::new(name).unwrap();
        unsafe { ffi::uno_get_solver_bool_option(self.ptr, name.as_ptr()) }
    }

    // logger

    pub fn set_logger_stream_callback(
        cb: UnoLoggerStreamCallback,
        user_data: *mut c_void,
    ) -> bool {
        unsafe { ffi::uno_set_logger_stream_callback(cb, user_data) }
    }

    pub fn reset_logger_stream() -> bool {
        unsafe { ffi::uno_reset_logger_stream() }
    }

    // callbacks

    pub fn set_callbacks(
        &self,
        notify_cb: Option<UnoNotifyAcceptableIterateCallback>,
        termination_cb: Option<UnoTerminationCallback>,
        user_data: *mut c_void,
    ) -> bool {
        unsafe { ffi::uno_set_solver_callbacks(self.ptr, notify_cb, termination_cb, user_data) }
    }

    // solve

    /// Optimizes `model` and returns a `Solution`.
    pub fn optimize(&self, model: &Model) -> Solution {
        unsafe { ffi::uno_optimize(self.ptr, model.ptr) }

        let _n = model.number_variables as usize;
        // We need to know the number of constraints to fill dual vectors.
        // Uno does not expose number_constraints directly from the solver,
        // but we can query per-component – we derive the count from the caller
        // by reading components until we get NaN (sentinel).  A simpler
        // approach: the caller provides the count.  For safety we just ask
        // the user to call `optimize_with_nc`.
        self.collect_solution(model.number_variables, 0)
    }

    /// Like `optimize` but also collects dual vectors.  Pass the number of
    /// constraints in `number_constraints`.
    pub fn optimize_with_nc(&self, model: &Model, number_constraints: usize) -> Solution {
        unsafe { ffi::uno_optimize(self.ptr, model.ptr) }
        self.collect_solution(model.number_variables, number_constraints as UnoInt)
    }

    fn collect_solution(&self, n: UnoInt, nc: UnoInt) -> Solution {
        let n_usize = n as usize;
        let nc_usize = nc as usize;

        let mut primals = vec![0.0f64; n_usize];
        let mut constraint_duals = vec![0.0f64; nc_usize];
        let mut lower_bound_duals = vec![0.0f64; n_usize];
        let mut upper_bound_duals = vec![0.0f64; n_usize];

        unsafe {
            ffi::uno_get_primal_solution(self.ptr, primals.as_mut_ptr());
            if nc_usize > 0 {
                ffi::uno_get_constraint_dual_solution(self.ptr, constraint_duals.as_mut_ptr());
            }
            ffi::uno_get_lower_bound_dual_solution(self.ptr, lower_bound_duals.as_mut_ptr());
            ffi::uno_get_upper_bound_dual_solution(self.ptr, upper_bound_duals.as_mut_ptr());
        }

        Solution {
            optimization_status: unsafe { ffi::uno_get_optimization_status(self.ptr) },
            solution_status: unsafe { ffi::uno_get_solution_status(self.ptr) },
            objective: unsafe { ffi::uno_get_solution_objective(self.ptr) },
            primals,
            constraint_duals,
            lower_bound_duals,
            upper_bound_duals,
            primal_feasibility: unsafe { ffi::uno_get_solution_primal_feasibility(self.ptr) },
            stationarity: unsafe { ffi::uno_get_solution_stationarity(self.ptr) },
            complementarity: unsafe { ffi::uno_get_solution_complementarity(self.ptr) },
            iterations: unsafe { ffi::uno_get_number_iterations(self.ptr) },
            cpu_time: unsafe { ffi::uno_get_cpu_time(self.ptr) },
        }
    }
}

impl Default for Solver {
    fn default() -> Self {
        Self::new()
    }
}
