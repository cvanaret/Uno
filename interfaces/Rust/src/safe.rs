// uno-rs: safe, idiomatic layer on top of the raw C callbacks.
//
// This module lets users provide ordinary Rust closures (or functions)
// instead of `unsafe extern "C"` callbacks. The crate stores the closures
// and installs generic trampolines that recover them from Uno's `user_data`
// pointer and call them with safe `&[f64]` / `&mut [f64]` slices.

use std::os::raw::{c_char, c_void};

use crate::ffi;
use crate::uno_int;

// ─────────────────────────────────────────────────────────────────────────
// Error type
// ─────────────────────────────────────────────────────────────────────────

/// Error returned by the safe builder API when an underlying Uno call fails.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnoError {
    /// Name of the underlying Uno C function that returned `false`.
    pub operation: &'static str,
}

impl std::fmt::Display for UnoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Uno operation `{}` failed", self.operation)
    }
}

impl std::error::Error for UnoError {}

type Result<T> = std::result::Result<T, UnoError>;

fn check(ok: bool, operation: &'static str) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(UnoError { operation })
    }
}

// ─────────────────────────────────────────────────────────────────────────
// User callback signatures (safe)
//
// Each callback returns `Result<(), ()>`: `Ok(())` maps to the Uno success
// code (0) and `Err(())` maps to a non-zero evaluation-error code, which is
// how the C API signals "could not evaluate at this point".
// ─────────────────────────────────────────────────────────────────────────

/// `f(x) -> objective_value`.
///
/// Inputs: `x` (length `n`). Output: write the scalar objective into `obj`.
pub type ObjectiveFn = dyn FnMut(&[f64], &mut f64) -> std::result::Result<(), ()>;

/// `∇f(x)`.
///
/// Inputs: `x` (length `n`). Output: fill `gradient` (length `n`).
pub type ObjectiveGradientFn = dyn FnMut(&[f64], &mut [f64]) -> std::result::Result<(), ()>;

/// `c(x)`.
///
/// Inputs: `x` (length `n`). Output: fill `constraints` (length `m`).
pub type ConstraintsFn = dyn FnMut(&[f64], &mut [f64]) -> std::result::Result<(), ()>;

/// Constraint Jacobian values in the COO order fixed by the sparsity pattern.
///
/// Inputs: `x` (length `n`). Output: fill `values` (length `nnz`).
pub type JacobianFn = dyn FnMut(&[f64], &mut [f64]) -> std::result::Result<(), ()>;

/// Lagrangian Hessian values (lower/upper triangle as declared).
///
/// Inputs: `x` (length `n`), `objective_multiplier` (ρ), `multipliers`
/// (length `m`). Output: fill `values` (length `nnz`).
pub type LagrangianHessianFn =
    dyn FnMut(&[f64], f64, &[f64], &mut [f64]) -> std::result::Result<(), ()>;

// ─────────────────────────────────────────────────────────────────────────
// Storage for the boxed closures.
//
// We keep every closure the user registers alive for as long as the model
// lives, in a single heap allocation whose address we hand to Uno as
// `user_data`. The trampolines below cast that pointer back to `&mut Callbacks`.
// ─────────────────────────────────────────────────────────────────────────

#[derive(Default)]
struct Callbacks {
    objective: Option<Box<ObjectiveFn>>,
    objective_gradient: Option<Box<ObjectiveGradientFn>>,
    constraints: Option<Box<ConstraintsFn>>,
    jacobian: Option<Box<JacobianFn>>,
    hessian: Option<Box<LagrangianHessianFn>>,
}

#[inline]
fn ret_code(r: std::result::Result<(), ()>) -> uno_int {
    match r {
        Ok(()) => 0,
        // Non-zero signals an evaluation error to Uno.
        Err(()) => 1,
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Trampolines: `extern "C"` functions matching the raw Uno signatures.
//
// SAFETY: Uno calls these with the `user_data` pointer we registered, which
// always points to the `Callbacks` box owned by the `Problem`. The box
// outlives every call because `optimize` borrows `&Problem`. The closures are
// only ever invoked from inside `uno_optimize`, single-threaded, so the
// `&mut` reborrow below is sound (no aliasing across calls).
// ─────────────────────────────────────────────────────────────────────────

unsafe fn callbacks<'a>(user_data: *mut c_void) -> &'a mut Callbacks {
    debug_assert!(!user_data.is_null(), "Uno called a trampoline with null user_data");
    &mut *(user_data as *mut Callbacks)
}

unsafe extern "C" fn objective_trampoline(
    n: uno_int,
    x: *const f64,
    obj: *mut f64,
    user_data: *mut c_void,
) -> uno_int {
    let cb = callbacks(user_data);
    let f = cb.objective.as_mut().expect("objective callback missing");
    let x = std::slice::from_raw_parts(x, n as usize);
    ret_code(f(x, &mut *obj))
}

unsafe extern "C" fn objective_gradient_trampoline(
    n: uno_int,
    x: *const f64,
    gradient: *mut f64,
    user_data: *mut c_void,
) -> uno_int {
    let cb = callbacks(user_data);
    let f = cb
        .objective_gradient
        .as_mut()
        .expect("objective-gradient callback missing");
    let x = std::slice::from_raw_parts(x, n as usize);
    let g = std::slice::from_raw_parts_mut(gradient, n as usize);
    ret_code(f(x, g))
}

unsafe extern "C" fn constraints_trampoline(
    n: uno_int,
    m: uno_int,
    x: *const f64,
    values: *mut f64,
    user_data: *mut c_void,
) -> uno_int {
    let cb = callbacks(user_data);
    let f = cb.constraints.as_mut().expect("constraints callback missing");
    let x = std::slice::from_raw_parts(x, n as usize);
    let c = std::slice::from_raw_parts_mut(values, m as usize);
    ret_code(f(x, c))
}

unsafe extern "C" fn jacobian_trampoline(
    n: uno_int,
    nnz: uno_int,
    x: *const f64,
    values: *mut f64,
    user_data: *mut c_void,
) -> uno_int {
    let cb = callbacks(user_data);
    let f = cb.jacobian.as_mut().expect("jacobian callback missing");
    let x = std::slice::from_raw_parts(x, n as usize);
    let v = std::slice::from_raw_parts_mut(values, nnz as usize);
    ret_code(f(x, v))
}

unsafe extern "C" fn hessian_trampoline(
    n: uno_int,
    m: uno_int,
    nnz: uno_int,
    x: *const f64,
    objective_multiplier: f64,
    multipliers: *const f64,
    values: *mut f64,
    user_data: *mut c_void,
) -> uno_int {
    let cb = callbacks(user_data);
    let f = cb.hessian.as_mut().expect("hessian callback missing");
    let x = std::slice::from_raw_parts(x, n as usize);
    let y = std::slice::from_raw_parts(multipliers, m as usize);
    let h = std::slice::from_raw_parts_mut(values, nnz as usize);
    ret_code(f(x, objective_multiplier, y, h))
}

// ─────────────────────────────────────────────────────────────────────────
// Problem: the safe model builder.
// ─────────────────────────────────────────────────────────────────────────

/// A safe, idiomatic optimization model.
///
/// Build it with [`Problem::new`], attach an objective / constraints / Hessian
/// using ordinary Rust closures, then hand it to [`crate::Solver::solve`].
///
/// The closures are stored on the heap inside the `Problem` and kept alive
/// until the `Problem` is dropped. The FFI plumbing (function pointers,
/// `user_data`, raw slices) is handled internally.
pub struct Problem {
    ptr: *mut c_void,
    number_variables: uno_int,
    number_constraints: uno_int,
    // Boxed twice: the outer Box gives a stable address to register as
    // `user_data`; the inner struct owns the user's closures. Heap-stable for
    // the lifetime of the Problem.
    callbacks: Box<Callbacks>,
}

impl Drop for Problem {
    fn drop(&mut self) {
        unsafe { ffi::uno_destroy_model(self.ptr) }
        // `self.callbacks` is freed here, after the model that referenced it.
    }
}

impl Problem {
    /// Creates a new model with per-variable bounds.
    ///
    /// Use `f64::NEG_INFINITY` / `f64::INFINITY` for unbounded variables.
    /// `problem_type` is one of `"NLP"`, `"QP"`, `"LP"`.
    pub fn new(
        problem_type: &str,
        lower_bounds: &[f64],
        upper_bounds: &[f64],
        base_indexing: uno_int,
    ) -> Result<Self> {
        assert_eq!(
            lower_bounds.len(),
            upper_bounds.len(),
            "lower and upper bound arrays must have the same length"
        );
        let n = lower_bounds.len() as uno_int;
        let pt = std::ffi::CString::new(problem_type).expect("problem_type contained a NUL byte");
        let ptr = unsafe {
            ffi::uno_create_model(
                pt.as_ptr(),
                n,
                lower_bounds.as_ptr(),
                upper_bounds.as_ptr(),
                base_indexing,
            )
        };
        if ptr.is_null() {
            return Err(UnoError { operation: "uno_create_model" });
        }

        let mut problem = Self {
            ptr,
            number_variables: n,
            number_constraints: 0,
            callbacks: Box::new(Callbacks::default()),
        };
        // Register the (stable) address of our callbacks box as Uno's user_data.
        let ud = problem.callbacks.as_mut() as *mut Callbacks as *mut c_void;
        check(
            unsafe { ffi::uno_set_user_data(problem.ptr, ud) },
            "uno_set_user_data",
        )?;
        Ok(problem)
    }

    /// Number of variables in the model.
    pub fn number_variables(&self) -> usize {
        self.number_variables as usize
    }

    /// Number of constraints registered so far.
    pub fn number_constraints(&self) -> usize {
        self.number_constraints as usize
    }

    /// Sets the objective from Rust closures.
    ///
    /// * `objective` receives `x` and writes the scalar value.
    /// * `gradient` receives `x` and fills the length-`n` gradient slice.
    ///
    /// Either closure may return `Err(())` to signal an evaluation error.
    pub fn set_objective<F, G>(
        &mut self,
        sense: uno_int,
        objective: F,
        gradient: G,
    ) -> Result<&mut Self>
    where
        F: FnMut(&[f64], &mut f64) -> std::result::Result<(), ()> + 'static,
        G: FnMut(&[f64], &mut [f64]) -> std::result::Result<(), ()> + 'static,
    {
        self.callbacks.objective = Some(Box::new(objective));
        self.callbacks.objective_gradient = Some(Box::new(gradient));
        check(
            unsafe {
                ffi::uno_set_objective(
                    self.ptr,
                    sense,
                    objective_trampoline,
                    objective_gradient_trampoline,
                )
            },
            "uno_set_objective",
        )?;
        Ok(self)
    }

    /// Sets the constraints and their Jacobian from Rust closures.
    ///
    /// The Jacobian sparsity is given in COO form by `jacobian_row_indices`
    /// and `jacobian_col_indices`; `jacobian` must fill its `values` slice in
    /// that same order.
    pub fn set_constraints<C, J>(
        &mut self,
        constraints: C,
        lower_bounds: &[f64],
        upper_bounds: &[f64],
        jacobian_row_indices: &[uno_int],
        jacobian_col_indices: &[uno_int],
        jacobian: J,
    ) -> Result<&mut Self>
    where
        C: FnMut(&[f64], &mut [f64]) -> std::result::Result<(), ()> + 'static,
        J: FnMut(&[f64], &mut [f64]) -> std::result::Result<(), ()> + 'static,
    {
        assert_eq!(
            lower_bounds.len(),
            upper_bounds.len(),
            "constraint lower/upper bound arrays must have the same length"
        );
        assert_eq!(
            jacobian_row_indices.len(),
            jacobian_col_indices.len(),
            "jacobian row/col index arrays must have the same length"
        );
        let nc = lower_bounds.len() as uno_int;
        let nnz = jacobian_row_indices.len() as uno_int;

        self.callbacks.constraints = Some(Box::new(constraints));
        self.callbacks.jacobian = Some(Box::new(jacobian));

        check(
            unsafe {
                ffi::uno_set_constraints(
                    self.ptr,
                    nc,
                    constraints_trampoline,
                    lower_bounds.as_ptr(),
                    upper_bounds.as_ptr(),
                    nnz,
                    jacobian_row_indices.as_ptr(),
                    jacobian_col_indices.as_ptr(),
                    jacobian_trampoline,
                )
            },
            "uno_set_constraints",
        )?;
        self.number_constraints = nc;
        Ok(self)
    }

    /// Sets the Lagrangian Hessian from a Rust closure.
    ///
    /// The sparsity is given in COO form by `row_indices` / `col_indices` over
    /// the declared triangle; the closure fills `values` in that order.
    /// Signature: `|x, objective_multiplier, multipliers, values|`.
    pub fn set_lagrangian_hessian<H>(
        &mut self,
        triangular_part: c_char,
        row_indices: &[uno_int],
        col_indices: &[uno_int],
        hessian: H,
    ) -> Result<&mut Self>
    where
        H: FnMut(&[f64], f64, &[f64], &mut [f64]) -> std::result::Result<(), ()> + 'static,
    {
        assert_eq!(
            row_indices.len(),
            col_indices.len(),
            "hessian row/col index arrays must have the same length"
        );
        let nnz = row_indices.len() as uno_int;

        self.callbacks.hessian = Some(Box::new(hessian));

        check(
            unsafe {
                ffi::uno_set_lagrangian_hessian(
                    self.ptr,
                    nnz,
                    triangular_part,
                    row_indices.as_ptr(),
                    col_indices.as_ptr(),
                    hessian_trampoline,
                )
            },
            "uno_set_lagrangian_hessian",
        )?;
        Ok(self)
    }

    /// Sets the Lagrangian sign convention.
    pub fn set_lagrangian_sign_convention(&mut self, convention: uno_int) -> Result<&mut Self> {
        check(
            unsafe { ffi::uno_set_lagrangian_sign_convention(self.ptr, convention) },
            "uno_set_lagrangian_sign_convention",
        )?;
        Ok(self)
    }

    /// Sets the initial primal iterate (length `n`).
    pub fn set_initial_primal_iterate(&mut self, x0: &[f64]) -> Result<&mut Self> {
        assert_eq!(
            x0.len(),
            self.number_variables as usize,
            "initial primal iterate has wrong length"
        );
        check(
            unsafe { ffi::uno_set_initial_primal_iterate(self.ptr, x0.as_ptr()) },
            "uno_set_initial_primal_iterate",
        )?;
        Ok(self)
    }

    /// Sets the initial dual iterate.
    pub fn set_initial_dual_iterate(&mut self, y0: &[f64]) -> Result<&mut Self> {
        check(
            unsafe { ffi::uno_set_initial_dual_iterate(self.ptr, y0.as_ptr()) },
            "uno_set_initial_dual_iterate",
        )?;
        Ok(self)
    }

    /// Raw model pointer, for use by the solver. Internal.
    pub(crate) fn as_ptr(&self) -> *mut c_void {
        self.ptr
    }
}
