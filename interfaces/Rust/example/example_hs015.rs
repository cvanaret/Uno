// Rust port of example_hs015.c
// HS015: min  100*(x1 - x0^2)^2 + (1-x0)^2
//        s.t. x0*x1 >= 1
//             x0 + x1^2 >= 0
//             x0 <= 0.5
// Solution: x = (0.5, 2.0)

use uno_rs::{
    get_version, uno_int, Model, Solver,
    UNO_FEASIBLE_KKT_POINT, UNO_LOWER_TRIANGLE, UNO_MINIMIZE,
    UNO_MULTIPLIER_NEGATIVE, UNO_SUCCESS, UNO_ZERO_BASED_INDEXING,
};
use std::os::raw::c_void;

// ── Objective: f(x) = 100*(x1 - x0^2)^2 + (1 - x0)^2 ──────────────────────

unsafe extern "C" fn objective_function(
    _n: uno_int, x: *const f64, obj: *mut f64, _data: *mut c_void,
) -> uno_int {
    let x = std::slice::from_raw_parts(x, 2);
    *obj = 100.0 * (x[1] - x[0] * x[0]).powi(2) + (1.0 - x[0]).powi(2);
    0
}

unsafe extern "C" fn objective_gradient(
    _n: uno_int, x: *const f64, g: *mut f64, _data: *mut c_void,
) -> uno_int {
    let x = std::slice::from_raw_parts(x, 2);
    let g = std::slice::from_raw_parts_mut(g, 2);
    g[0] = 400.0 * x[0].powi(3) - 400.0 * x[0] * x[1] + 2.0 * x[0] - 2.0;
    g[1] = 200.0 * (x[1] - x[0] * x[0]);
    0
}

// ── Constraints: c0 = x0*x1,  c1 = x0 + x1^2 ───────────────────────────────

unsafe extern "C" fn constraint_functions(
    _n: uno_int, _nc: uno_int, x: *const f64, cv: *mut f64, _data: *mut c_void,
) -> uno_int {
    let x = std::slice::from_raw_parts(x, 2);
    let cv = std::slice::from_raw_parts_mut(cv, 2);
    cv[0] = x[0] * x[1];
    cv[1] = x[0] + x[1] * x[1];
    0
}

// ── Jacobian (COO, 4 nonzeros) ───────────────────────────────────────────────
// Row 0 (c0): ∂/∂x0 = x1  (col 0),  ∂/∂x1 = x0   (col 1)
// Row 1 (c1): ∂/∂x0 = 1   (col 0),  ∂/∂x1 = 2*x1 (col 1)

unsafe extern "C" fn jacobian_values(
    _n: uno_int, _nnz: uno_int, x: *const f64, jv: *mut f64, _data: *mut c_void,
) -> uno_int {
    let x = std::slice::from_raw_parts(x, 2);
    let jv = std::slice::from_raw_parts_mut(jv, 4);
    jv[0] = x[1];        // (row 0, col 0)
    jv[1] = 1.0;         // (row 1, col 0)
    jv[2] = x[0];        // (row 0, col 1)
    jv[3] = 2.0 * x[1];  // (row 1, col 1)
    0
}

// ── Lagrangian Hessian – lower triangle, 3 nonzeros ─────────────────────────
// L = rho*f(x) - y^T c(x)   (UNO_MULTIPLIER_NEGATIVE convention)
// H00 = rho*(1200*x0^2 - 400*x1 + 2)
// H10 = -400*rho*x0 - y0
// H11 = 200*rho - 2*y1

unsafe extern "C" fn lagrangian_hessian(
    _n: uno_int, _nc: uno_int, _nnz: uno_int,
    x: *const f64, rho: f64, y: *const f64, hv: *mut f64,
    _data: *mut c_void,
) -> uno_int {
    let x = std::slice::from_raw_parts(x, 2);
    let y = std::slice::from_raw_parts(y, 2);
    let hv = std::slice::from_raw_parts_mut(hv, 3);
    hv[0] = rho * (1200.0 * x[0] * x[0] - 400.0 * x[1] + 2.0);
    hv[1] = -400.0 * rho * x[0] - y[0];
    hv[2] = 200.0 * rho - 2.0 * y[1];
    0
}

// ────────────────────────────────────────────────────────────────────────────

fn print_vec(label: &str, v: &[f64]) {
    print!("{}: ", label);
    for val in v { print!("{:.6} ", val); }
    println!();
}

fn main() {
    let (maj, min, pat) = get_version();
    println!("Uno v{}.{}.{}", maj, min, pat);

    // ── Build model ────────────────────────────────────────────────────────────
    let var_lb = [f64::NEG_INFINITY, f64::NEG_INFINITY];
    let var_ub = [0.5_f64, f64::INFINITY];
    let mut model = Model::new("NLP", &var_lb, &var_ub, UNO_ZERO_BASED_INDEXING);

    assert!(model.set_objective(UNO_MINIMIZE, objective_function, objective_gradient));

    let con_lb = [1.0_f64, 0.0_f64];
    let con_ub = [f64::INFINITY, f64::INFINITY];
    // COO sparsity pattern: rows then cols, column-major ordering
    let jac_rows: [uno_int; 4] = [0, 1, 0, 1];
    let jac_cols: [uno_int; 4] = [0, 0, 1, 1];
    assert!(model.set_constraints(
        constraint_functions,
        &con_lb, &con_ub,
        &jac_rows, &jac_cols,
        jacobian_values,
    ));

    assert!(model.set_initial_primal_iterate(&[-2.0_f64, 1.0_f64]));

    // ── Build solver ───────────────────────────────────────────────────────────
    let solver = Solver::new();
    // "ipopt" preset: interior-point method — only needs MUMPS (no QP solver)
    solver.set_preset("ipopt");
    solver.set_bool_option("print_solution", true);

    // ── Run 1: L-BFGS Hessian (default for NLP when no Hessian provided) ──────
    println!("\n=== Run 1: L-BFGS Hessian ===");
    let sol1 = solver.optimize_with_nc(&model, 2);
    println!("Optimization status : {}", sol1.optimization_status_str());
    println!("Solution status     : {}", sol1.solution_status_str());
    println!("Objective           : {:.10}", sol1.objective);
    assert_eq!(sol1.optimization_status, UNO_SUCCESS,
               "Run 1 optimization failed");
    assert_eq!(sol1.solution_status, UNO_FEASIBLE_KKT_POINT,
               "Run 1 did not find a KKT point");

    // ── Run 2: Exact Hessian ───────────────────────────────────────────────────
    println!("\n=== Run 2: Exact Hessian ===");
    let hess_rows: [uno_int; 3] = [0, 1, 1];
    let hess_cols: [uno_int; 3] = [0, 0, 1];
    assert!(model.set_lagrangian_hessian(
        UNO_LOWER_TRIANGLE,
        &hess_rows, &hess_cols,
        lagrangian_hessian,
    ));
    assert!(model.set_lagrangian_sign_convention(UNO_MULTIPLIER_NEGATIVE));
    solver.set_string_option("hessian_model", "exact");

    let sol2 = solver.optimize_with_nc(&model, 2);
    println!("Optimization status : {}", sol2.optimization_status_str());
    println!("Solution status     : {}", sol2.solution_status_str());
    assert_eq!(sol2.optimization_status, UNO_SUCCESS,
               "Run 2 optimization failed");
    assert_eq!(sol2.solution_status, UNO_FEASIBLE_KKT_POINT,
               "Run 2 did not find a KKT point");

    println!("Objective           : {:.10}", sol2.objective);
    print_vec("Primal solution    ", &sol2.primals);
    print_vec("Constraint duals   ", &sol2.constraint_duals);
    print_vec("Lower-bound duals  ", &sol2.lower_bound_duals);
    print_vec("Upper-bound duals  ", &sol2.upper_bound_duals);
    println!("Primal feasibility  : {:.2e}", sol2.primal_feasibility);
    println!("Stationarity        : {:.2e}", sol2.stationarity);
    println!("Complementarity     : {:.2e}", sol2.complementarity);
    println!("Iterations          : {}", sol2.iterations);
    println!("CPU time            : {:.4}s", sol2.cpu_time);

    // ── Verify x ≈ (0.5, 2.0) ─────────────────────────────────────────────────
    let tol = 1e-4;
    let x = &sol2.primals;
    println!("\n=== Verification ===");
    println!("Expected : x0 = 0.5,  x1 = 2.0");
    println!("Got      : x0 = {:.8}, x1 = {:.8}", x[0], x[1]);
    assert!((x[0] - 0.5).abs() < tol, "x[0] = {} not near 0.5", x[0]);
    assert!((x[1] - 2.0).abs() < tol, "x[1] = {} not near 2.0", x[1]);
    println!("✓  Solution x = (0.5, 2.0) confirmed within tolerance {}", tol);
}
