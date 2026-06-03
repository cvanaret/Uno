## Uno's Rust interface: how to use Uno from Rust

Uno's Rust interface (`uno_rs`) allows you to solve an optimization model described by callback functions.

There are two ways to use the crate:

- **The idiomatic API (recommended)** lets you describe your problem with ordinary Rust closures or functions. The crate handles the FFI layer (function pointers, `user_data`, raw slices) for you. See [example_hs015_idiomatic.rs](example/example_hs015_idiomatic.rs).
- **The raw API** mirrors the C interface one-to-one: you pass `unsafe extern "C"` callbacks and raw pointers yourself. See [example_hs015.rs](example/example_hs015.rs). Use this only if you need full control over the FFI boundary.

Both APIs share the same `Solver`, options, and `Solution` types.

### Idiomatic API (closures)

Import the safe builder:

```rust
use uno_rs::{uno_int, Problem, Solver, UNO_MINIMIZE, UNO_ZERO_BASED_INDEXING};
```

Build a [`Problem`] and attach an objective, constraints and (optionally) a Hessian using closures. Each closure takes safe `&[f64]` inputs and writes into a `&mut [f64]` (or `&mut f64`) output, returning `Ok(())` on success or `Err(())` to signal an evaluation error at the current point:

```rust
let var_lb = [f64::NEG_INFINITY, f64::NEG_INFINITY];
let var_ub = [0.5, f64::INFINITY];
let mut problem = Problem::new("NLP", &var_lb, &var_ub, UNO_ZERO_BASED_INDEXING)?;

// Objective and its gradient.
problem.set_objective(
    UNO_MINIMIZE,
    |x, obj| { *obj = 100.0 * (x[1] - x[0]*x[0]).powi(2) + (1.0 - x[0]).powi(2); Ok(()) },
    |x, g| {
        g[0] = 400.0*x[0].powi(3) - 400.0*x[0]*x[1] + 2.0*x[0] - 2.0;
        g[1] = 200.0 * (x[1] - x[0]*x[0]);
        Ok(())
    },
)?;

// Constraints c0 = x0*x1, c1 = x0 + x1^2, plus their Jacobian (COO order).
let con_lb = [1.0, 0.0];
let con_ub = [f64::INFINITY, f64::INFINITY];
let jac_rows = [0, 1, 0, 1];
let jac_cols = [0, 0, 1, 1];
problem.set_constraints(
    |x, c| { c[0] = x[0]*x[1]; c[1] = x[0] + x[1]*x[1]; Ok(()) },
    &con_lb, &con_ub, &jac_rows, &jac_cols,
    |x, jv| { jv[0] = x[1]; jv[1] = 1.0; jv[2] = x[0]; jv[3] = 2.0*x[1]; Ok(()) },
)?;

problem.set_initial_primal_iterate(&[-2.0, 1.0])?;

let solver = Solver::new();
solver.set_preset("ipopt");
let solution = solver.solve(&problem);   // constraint count is taken from the Problem
```

Because these are real closures, they may capture state from their environment — something the raw `extern "C"` function-pointer callbacks cannot do.

The builder methods return `Result<&mut Self, UnoError>`, so failures surface as ordinary Rust errors (and the `&mut Self` return makes chaining with `?` convenient). `Solver::solve` takes a `&Problem` and reads the constraint count from it, so dual vectors always come back with the right length.

The rest of this document describes the **raw API**.

Start by importing from the crate:

```rust
use uno_rs::{uno_int, Model, Solver, UNO_MINIMIZE, UNO_ZERO_BASED_INDEXING};
use std::os::raw::c_void;
```

### Building an optimization model

Building an optimization model is incremental and starts with the information about the variables:

```rust
let model = Model::new(problem_type, &variables_lower_bounds, &variables_upper_bounds, base_indexing);
```

The following optional elements can be added to the model separately:
- the objective function (and its gradient). It is 0 otherwise;
```rust
model.set_objective(optimization_sense, objective_function, objective_gradient);
```
- constraint functions (and their Jacobian);
```rust
model.set_constraints(constraint_functions, &constraints_lower_bounds, &constraints_upper_bounds,
   &jacobian_row_indices, &jacobian_column_indices, jacobian);
```
- the Lagrangian Hessian;
```rust
model.set_lagrangian_hessian(hessian_triangular_part, &hessian_row_indices, &hessian_column_indices, lagrangian_hessian);
```
- a Jacobian operator (performs Jacobian-vector products);
```rust
model.set_jacobian_operator(jacobian_operator);
```
- a Jacobian-transposed operator (performs Jacobian-transposed-vector products);
```rust
model.set_jacobian_transposed_operator(jacobian_transposed_operator);
```
- a Hessian operator (performs Hessian-vector products);
```rust
model.set_lagrangian_hessian_operator(lagrangian_hessian_operator);
```
- a Lagrangian sign convention (default is `UNO_MULTIPLIER_NEGATIVE`);
```rust
model.set_lagrangian_sign_convention(lagrangian_sign_convention);
```
- user data of an arbitrary type (`*mut c_void`);
```rust
model.set_user_data(user_data);
```
- an initial primal point;
```rust
model.set_initial_primal_iterate(&initial_primal_iterate);
```
- an initial dual point.
```rust
model.set_initial_dual_iterate(&initial_dual_iterate);
```

*Each of these functions returns a `bool` that is `true` upon success and `false` upon failure.*

The memory for the model is allocated by Uno and freed automatically when the `Model` is dropped.

### Creating an instance of the Uno solver

Create an instance of the Uno solver with:
```rust
let solver = Solver::new();
```

The memory for the solver is allocated by Uno and freed automatically when the `Solver` is dropped.

### Passing options to the Uno solver

Options can be passed to the Uno solver:
```rust
solver.set_integer_option("max_iterations", 1000);
solver.set_double_option("primal_tolerance", 1.0e-6);
solver.set_bool_option("print_solution", true);
solver.set_string_option("hessian_model", "exact");
```

Loading options from a file (overwrites existing options):
```rust
solver.load_option_file("uno.opt");
```

Getting typed value of an option from the Uno solver:
```rust
solver.get_string_option("hessian_model");
solver.get_bool_option("print_solution");
```

Setting a preset has Uno mimic an existing solver:
```rust
solver.set_preset("filtersqp");
```

### Setting solver callbacks

Setting the user callbacks to the Uno solver:
```rust
solver.set_callbacks(notify_acceptable_iterate_callback, user_termination_callback, user_data);
```

Setting the logger stream callback:
```rust
Solver::set_logger_stream_callback(logger_stream_callback, user_data);
```

and reset the logger stream to the standard output:
```rust
Solver::reset_logger_stream();
```

### Solving the model

The model can then be solved by Uno:
```rust
let solution = solver.optimize(&model);
// or, to also retrieve dual vectors:
let solution = solver.optimize_with_nc(&model, number_constraints);
```

### Inspecting the result

A `Solution` struct allows you to inspect the result of the optimization:
- the optimization status (`UNO_SUCCESS`, `UNO_ITERATION_LIMIT`, `UNO_TIME_LIMIT`, `UNO_EVALUATION_ERROR`, `UNO_ALGORITHMIC_ERROR`): `solution.optimization_status`
- the solution status (`UNO_NOT_OPTIMAL`, `UNO_FEASIBLE_KKT_POINT`, `UNO_FEASIBLE_FJ_POINT`, `UNO_INFEASIBLE_STATIONARY_POINT`, `UNO_FEASIBLE_SMALL_STEP`, `UNO_INFEASIBLE_SMALL_STEP`, `UNO_UNBOUNDED`): `solution.solution_status`
- the objective value of the solution: `solution.objective`
- the primal solution: `solution.primals`
- the dual solution associated with the general constraints: `solution.constraint_duals`
- the dual solution associated with the lower bounds: `solution.lower_bound_duals`
- the dual solution associated with the upper bounds: `solution.upper_bound_duals`
- the primal feasibility residual at the solution: `solution.primal_feasibility`
- the stationarity residual at the solution: `solution.stationarity`
- the complementarity residual at the solution: `solution.complementarity`
- the number of (outer) iterations: `solution.iterations`
- the CPU time (in seconds): `solution.cpu_time`
