%% HS015 problem
clc, clear

%% Optimization model
% Problem type: 'L' = Linear, 'Q' = Quadratic, 'N' = Nonlinear
model.problem_type = 'N';
% Vector indexing
model.base_indexing = 1;
% Number of variables
model.number_variables = 2;
% Variable bounds
model.variables_lower_bounds = [-inf, -inf]; 
model.variables_upper_bounds = [0.5, inf];

% Optimization sense: 1 = Minimize, -1 = Maximize
model.optimization_sense = 1;
% Objective function handle: objective = objective_function(x)
model.objective_function = @objective_function; 
% Objective gradient handle: gradient = objective_gradient(x)
model.objective_gradient = @objective_gradient;

% Number of constraints
model.number_constraints = 2;
% Constraint bounds
model.constraints_lower_bounds = [1.0, 0.0];
model.constraints_upper_bounds = [inf, inf];
% Constraint function handle: constraints = constrain_function(x)
model.constraint_function = @constraint_function;
% Constraint jacobian handle: jacobian = constraint_jacobian(x)
model.constraint_jacobian = @constraint_jacobian;
% Constraint sparsity pattern (base_indexing-based)
model.number_jacobian_nonzeros = 4;
model.jacobian_row_indices = [1 2 1 2];
model.jacobian_column_indices = [1 1 2 2];

% Lagrangian sign convention: 1 = rho*f(x) + y^T c(x), -1 = rho*f(x) - y^T c(x)
model.lagrangian_sign_convention = -1;
% Hessian triangular part: 'L' = lower, 'U' = upper
model.hessian_triangular_part = 'L';
% Lagrangian Hessian handle: hessian = lagrangian_hessian(x,rho,y)
model.lagrangian_hessian = @lagrangian_hessian;
% Lagrangian Hessian sparsity pattern (base_indexing-based)
model.number_hessian_nonzeros = 3;
model.hessian_row_indices = [1 2 2];
model.hessian_column_indices = [1 1 2];

% Jacobian operator handle: result = jacobian_operator(x,v)
% model.jacobian_operator = @jacobian_operator;

% Jacobian tranposed operator handle: result = jacobian_transposed_operator(x,v)
% model.jacobian_transposed_operator = @jacobian_transposed_operator;

% Hessian operator handle: result = lagrangian_hessian_operator(x,rho,y,v)
% model.lagrangian_hessian_operator = @lagrangian_hessian_operator;

% Initial primal point
model.initial_primal_iterate = [-2.0, 1.0];
% Dual primal point
model.initial_dual_iterate = [0, 0];

%% UNO solver options
% Preset can be possibly given: options = uno_options([preset])
options = uno_options(); % default preset
% options = uno_options('ipopt');

%% UNO solver callbacks
callbacks = struct();
% Logger stream handle: logger_stream_callback(str)
% callbacks.logger_stream_callback = @logger_stream_callback;
% Notify acceptable iterate handle: notify_acceptable_iterate_callback(x, yl, yb, y, rho, feas, stat, compl)
% callbacks.notify_acceptable_iterate_callback = @notify_acceptable_iterate_callback;
% Notify new primal handle: notify_new_primals(x)
% callbacks.notify_new_primals_callback = @notify_new_primals_callback;
% Notify new dual handle: notify_new_multipliers(x, yl, yb, y)
% callbacks.notify_new_multipliers_callback = @notify_new_multipliers_callback;
% User termination handle: terminate = user_termination_callback(x, yl, yb, y, rho, feas, stat, compl)
% callbacks.user_termination_callback = @user_termination_callback;

%% Solving the model
result = uno_optimize(model, options, callbacks);
% Clear mex functions
clear uno_optimize uno_options

%% Display the result
disp(result)

%% Model functions
% Objective function
function fval = objective_function(x)
    fval = 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
end

% Constraint function
function c = constraint_function(x)
    c = zeros(2,1);
    c(1) = x(1)*x(2);
    c(2) = x(1) + x(2)^2;
end

% Objective gradient
function grad = objective_gradient(x)
    grad = zeros(2,1);
    grad(1) = 400*x(1)^3 - 400*x(1)*x(2) + 2*x(1) - 2;
    grad(2) = 200*(x(2) - x(1)^2);
end

% Constraint jacobian
function J = constraint_jacobian(x)
    J = zeros(4,1);
    J(1) = x(2);
    J(2) = 1;
    J(3) = x(1);
    J(4) = 2*x(2);
end

% Lagrangian Hessian - lower triangular part
function H = lagrangian_hessian(x, objective_multiplier, multipliers)
    H = zeros(3,1);
    H(1) = objective_multiplier*(1200*x(1)^2 - 400*x(2) + 2);
    H(2) = -400*objective_multiplier*x(1) - multipliers(1);
    H(3) = 200*objective_multiplier - 2*multipliers(2);
end

% Lagrangian Hessian operator
function result = lagrangian_hessian_operator(x, objective_multiplier, multipliers, vector)
    h00 = objective_multiplier*(1200*x(1)^2 - 400*x(2) + 2);
    h10 = -400*objective_multiplier*x(1) - multipliers(1);
    h11 = 200*objective_multiplier - 2*multipliers(2);

    result = zeros(2,1);
    result(1) = h00*vector(1) + h10*vector(2);
    result(2) = h10*vector(1) + h11*vector(2);
end