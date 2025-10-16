## Uno's MATLAB interface: how to use Uno from MATLAB

Uno's MATLAB interface allows you to solve an optimization model described in MATLAB.
The MATLAB MEX functions `uno_optimize` and `uno_options` are compiled via the command `make unomex`.
An example is available in the file [example_hs015.m](example/example_hs015.m).

### Building an optimization model

Create optimization model as a MATLAB struct:

```matlab
% Problem type: 'L' = Linear, 'Q' = Quadratic, 'N' = Nonlinear
model.problem_type = ...;
% Vector indexing: 1 = Matlab-style, 0 = C-style
model.base_indexing = ...;
% Number of variables
model.number_variables = ...;
% Variable bounds
model.variables_lower_bounds = [ ... ]; 
model.variables_upper_bounds = [ ... ];
```

The following fields has also to be set:

- objective and its full gradient;
```matlab
% Optimization sense: 1 = Minimize, -1 = Maximize
model.optimization_sense = 1;
% Objective function handle: objective = objective_function(x)
model.objective_function = @objective_function; 
% Objective gradient handle: gradient = objective_gradient(x)
model.objective_gradient = @objective_gradient;
```

- constraint function and bounds and its sparse Jacobian matrix;
```matlab
% Number of constraints
model.number_constraints = ...;
% Constraint bounds
model.constraints_lower_bounds = [ ... ];
model.constraints_upper_bounds = [ ... ];
% Constraint function handle: constraints = constrain_function(x)
model.constraint_function = @constraint_function;
% Constraint jacobian handle: jacobian_values = constraint_jacobian(x)
model.constraint_jacobian = @constraint_jacobian;
% Constraint sparsity pattern (base_indexing-based)
model.number_jacobian_nonzeros = ...;
model.jacobian_row_indices = [ ... ];
model.jacobian_column_indices = [ ... ];
```

- sparse Lagrangian Hessian matrix;
```matlab
% Lagrangian sign convention: 1 = rho*f(x) + y^T c(x), -1 = rho*f(x) - y^T c(x)
model.lagrangian_sign_convention = ...;
% Hessian triangular part: 'L' = lower, 'U' = upper
model.hessian_triangular_part = 'L';
% Lagrangian Hessian handle: hessian_values = lagrangian_hessian(x,rho,y)
model.lagrangian_hessian = @lagrangian_hessian;
% Lagrangian Hessian sparsity pattern (base_indexing-based)
model.number_hessian_nonzeros = ...;
model.hessian_row_indices = [ ... ];
model.hessian_column_indices = [ ... ];
```

- a Jacobian operator (performs Jacobian-vector products);
```matlab
% Jacobian operator handle: result = jacobian_operator(x,v)
model.jacobian_operator = @jacobian_operator;
```

- a Jacobian-transposed operator (performs Jacobian-transposed-vector products);
```matlab
% Jacobian tranposed operator handle: result = jacobian_transposed_operator(x,v)
model.jacobian_transposed_operator = @jacobian_transposed_operator;
```

- a Hessian operator (performs Hessian-vector products);
```matlab
% Hessian operator handle: result = lagrangian_hessian_operator(x,rho,y,v)
model.lagrangian_hessian_operator = @lagrangian_hessian_operator;
```

- an initial primal and dual point;
```matlab
% Initial primal point
model.initial_primal_iterate = [ ... ];
% Dual primal point
model.initial_dual_iterate = [ ... ];
```

User data can be employed in the optimization model functions using function handle, for example
```matlab
model.objective_function = @(x) objective_function(x, user_data); 
```

### Defining UNO solver options

Optionally obtain the default options for the UNO solver as a MATLAB struct:
```matlab
% Preset can be possibly given: options = uno_options([preset])
options = uno_options();
```

Default preset is used if empty or no arguments are given

### Defining UNO solver callbacks

Optionally define UNO solver callbacks in a MATLAB struct:

- Logger stream callback (called when printing string to the stream)
```matlab
% Logger stream handle: logger_stream_callback(str)
callbacks.logger_stream_callback = @logger_stream_callback;
```

- Notify callback (called when acceptable iterate found)
```matlab
% Notify acceptable iterate handle: notify_acceptable_iterate_callback(x, yl, yb, y, rho, feas, stat, compl)
callbacks.notify_acceptable_iterate_callback = @notify_acceptable_iterate_callback;
```

- User termination callback (called at the major iteration)
```matlab
% User termination handle: terminate = user_termination_callback(x, yl, yb, y, rho, feas, stat, compl)
callbacks.user_termination_callback = @user_termination_callback;
```

### Solving the model

The model can then be solved by Uno, with optional input arguments `options` and `callbacks`:
```matlab
% Call UNO optimizer
% result = uno_optimize(model);
% result = uno_optimize(model, options);
result = uno_optimize(model, options, callbacks);
```

### Inspecting the result

To inspect the result of the optimization, read the fields of the `result` struct:
- the optimization status (`UNO_SUCCESS = 0`, `UNO_ITERATION_LIMIT = 1`, `UNO_TIME_LIMIT = 2`, `UNO_EVALUATION_ERROR = 3`, `UNO_ALGORITHMIC_ERROR = 4`, `UNO_USER_TERMINATION = 5`): `result.optimization_status`
- the solution status (`UNO_NOT_OPTIMAL = 0`, `UNO_FEASIBLE_KKT_POINT = 1`, `UNO_FEASIBLE_FJ_POINT = 2`, `UNO_INFEASIBLE_STATIONARY_POINT = 3`, `UNO_FEASIBLE_SMALL_STEP = 4`, `UNO_INFEASIBLE_SMALL_STEP = 5`, `UNO_UNBOUNDED = 6`): `result.solution_status`
- the objective value of the solution: `result.solution_objective`
- the primal solution: `result.primal_solution`
- the dual solution associated with the general constraints: `result.constraint_dual_solution`
- the dual solution associated with the lower bounds: `result.lower_bound_dual_solution`
- the dual solution associated with the upper bounds: `result.upper_bound_dual_solution`
- the primal feasibility residual at the solution: `result.solution_primal_feasibility`
- the stationarity residual at the solution: `result.solution_stationarity`
- the complementarity residual at the solution: `result.solution_complementarity`
- the number of (outer) iterations: `result.number_iterations`
- the CPU time (in seconds): `result.cpu_time`
- the number of objective evaluations: `result.number_objective_evaluations`
- the number of constraint evaluations: `result.number_constraint_evaluations`
- the number of objective gradient evaluations: `result.number_objective_gradient_evaluations`
- the number of Jacobian evaluations: `result.number_jacobian_evaluations`
- the number of Hessian evaluations: `result.number_hessian_evaluations`
- the number of subproblems solved: `result.number_subproblems_solved`