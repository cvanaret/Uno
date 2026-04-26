## Uno's Rust interface: how to use Uno from Rust

Uno's Rust interface (`uno_rs`) allows you to solve an optimization model described by callback functions.
An example is available in the file [example_hs015.rs](example/example_hs015.rs).

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
