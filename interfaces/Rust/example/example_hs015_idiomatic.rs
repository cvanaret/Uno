// Rust port of example_hs015.c, using the idiomatic `Problem` API.
// HS015: min  100*(x1 - x0^2)^2 + (1-x0)^2
//        s.t. x0*x1 >= 1
//             x0 + x1^2 >= 0
//             x0 <= 0.5
// Solution: x = (0.5, 2.0)
//
// Note how the objective, gradient, constraints, Jacobian and Hessian are
// ordinary Rust closures: no `unsafe extern "C"`, no raw pointers, no manual
// `from_raw_parts`. The crate installs the FFI trampolines internally.

use uno_rs::{
    get_version, uno_int, Problem, Solver,
    UNO_FEASIBLE_KKT_POINT, UNO_LOWER_TRIANGLE, UNO_MINIMIZE,
    UNO_MULTIPLIER_NEGATIVE, UNO_SUCCESS, UNO_ZERO_BASED_INDEXING,
};

fn print_vec(label: &str, v: &[f64]) {
    print!("{}: ", label);
    for val in v {
        print!("{:.6} ", val);
    }
    println!();
}

fn main() {
    let (maj, min, pat) = get_version();
    println!("Uno v{}.{}.{}", maj, min, pat);

    // ── Build model ──────────────────────────────────────────────────────────
    let var_lb = [f64::NEG_INFINITY, f64::NEG_INFINITY];
    let var_ub = [0.5_f64, f64::INFINITY];
    let mut problem =
        Problem::new("NLP", &var_lb, &var_ub, UNO_ZERO_BASED_INDEXING).expect("create model");

    // Objective f(x) = 100*(x1 - x0^2)^2 + (1 - x0)^2, plus its gradient.
    problem
        .set_objective(
            UNO_MINIMIZE,
            |x, obj| {
                *obj = 100.0 * (x[1] - x[0] * x[0]).powi(2) + (1.0 - x[0]).powi(2);
                Ok(())
            },
            |x, g| {
                g[0] = 400.0 * x[0].powi(3) - 400.0 * x[0] * x[1] + 2.0 * x[0] - 2.0;
                g[1] = 200.0 * (x[1] - x[0] * x[0]);
                Ok(())
            },
        )
        .expect("set objective");

    // Constraints c0 = x0*x1, c1 = x0 + x1^2, plus their Jacobian (COO).
    let con_lb = [1.0_f64, 0.0_f64];
    let con_ub = [f64::INFINITY, f64::INFINITY];
    let jac_rows: [uno_int; 4] = [0, 1, 0, 1];
    let jac_cols: [uno_int; 4] = [0, 0, 1, 1];
    problem
        .set_constraints(
            |x, c| {
                c[0] = x[0] * x[1];
                c[1] = x[0] + x[1] * x[1];
                Ok(())
            },
            &con_lb,
            &con_ub,
            &jac_rows,
            &jac_cols,
            |x, jv| {
                jv[0] = x[1]; // (row 0, col 0)
                jv[1] = 1.0; // (row 1, col 0)
                jv[2] = x[0]; // (row 0, col 1)
                jv[3] = 2.0 * x[1]; // (row 1, col 1)
                Ok(())
            },
        )
        .expect("set constraints");

    // Exact Lagrangian Hessian (lower triangle, 3 nonzeros).
    // L = rho*f(x) - y^T c(x)  (UNO_MULTIPLIER_NEGATIVE convention)
    let hess_rows: [uno_int; 3] = [0, 1, 1];
    let hess_cols: [uno_int; 3] = [0, 0, 1];
    problem
        .set_lagrangian_hessian(UNO_LOWER_TRIANGLE, &hess_rows, &hess_cols, |x, rho, y, hv| {
            hv[0] = rho * (1200.0 * x[0] * x[0] - 400.0 * x[1] + 2.0);
            hv[1] = -400.0 * rho * x[0] - y[0];
            hv[2] = 200.0 * rho - 2.0 * y[1];
            Ok(())
        })
        .expect("set hessian");
    problem
        .set_lagrangian_sign_convention(UNO_MULTIPLIER_NEGATIVE)
        .expect("set sign convention");

    problem
        .set_initial_primal_iterate(&[-2.0_f64, 1.0_f64])
        .expect("set initial iterate");

    // ── Build solver ───────────────────────────────────────────────────────────
    let solver = Solver::new();
    solver.set_preset("ipopt");
    solver.set_bool_option("print_solution", true);
    solver.set_string_option("hessian_model", "exact");

    // ── Solve ─────────────────────────────────────────────────────────────────
    let sol = solver.solve(&problem);
    println!("Optimization status : {}", sol.optimization_status_str());
    println!("Solution status     : {}", sol.solution_status_str());
    println!("Objective           : {:.10}", sol.objective);
    print_vec("Primal solution    ", &sol.primals);
    print_vec("Constraint duals   ", &sol.constraint_duals);

    assert_eq!(sol.optimization_status, UNO_SUCCESS, "optimization failed");
    assert_eq!(
        sol.solution_status, UNO_FEASIBLE_KKT_POINT,
        "did not find a KKT point"
    );

    println!("\n✓  Solved with the idiomatic closure-based API");
}
