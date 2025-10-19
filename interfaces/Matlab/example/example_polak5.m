%% Polak5 prolbem
clc; clear;

%% Optimization model
% Problem type: 'L' = Linear, 'Q' = Quadratic, 'N' = Nonlinear
model.problem_type = 'N';  
% Vector indexing
model.base_indexing = 1;   
% Number of variables
model.number_variables = 3;
% Variable bounds
model.variables_lower_bounds = [-inf, -inf, -inf];
model.variables_upper_bounds = [ inf,  inf,  inf];

% Optimization sense: 1 = Minimize, -1 = Maximize
model.optimization_sense = 1;
% Objective function handle: objective = objective_function(x)
model.objective_function = @objective_function;
% Objective gradient handle: gradient = objective_gradient(x)
model.objective_gradient = @objective_gradient;

% Number of constraints
model.number_constraints = 2;
% Constraint bounds
model.constraints_lower_bounds = [-inf, -inf];
model.constraints_upper_bounds = [0, 0];
% Constraint function handle: constraints = constrain_function(x)
model.constraint_function = @constraint_function;
% Constraint jacobian handle: jacobian = constraint_jacobian(x)
model.constraint_jacobian = @constraint_jacobian;
% Constraint sparsity pattern (base_indexing-based)
model.number_jacobian_nonzeros = 6;
model.jacobian_row_indices    = [1 1 1 2 2 2];
model.jacobian_column_indices = [1 2 3 1 2 3];

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
model.initial_primal_iterate = [0.1, 0.1, 0.0];
% Dual primal point
model.initial_dual_iterate   = [0, 0];

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

% User termination handle: terminate = user_termination_callback(x, yl, yb, y, rho, feas, stat, compl)
% callbacks.user_termination_callback = @user_termination_callback;

%% Solving the model
result = uno_optimize(model, options, callbacks);
% Clear mex functions
clear uno_optimize uno_options

%% Display the result
% disp(result)

%% Model functions
% Objective
function fval = objective_function(x)
    fval = x(3);
end

% Objective gradient
function grad = objective_gradient(~)
    grad = [0; 0; 1];
end

% Constraint functions
function c = constraint_function(x)
    c = zeros(2,1);
    c(1) = -x(3) + 3*x(1)^2 + 50*(x(1) - x(2)^4 - 1)^2;
    c(2) = -x(3) + 3*x(1)^2 + 50*(x(1) - x(2)^4 + 1)^2;
end

% Constraint Jacobian
function J = constraint_jacobian(x)
    J = zeros(6,1);
    J(1) = 6*x(1) + 100*(x(1) - x(2)^4 - 1);
    J(2) = -400*x(2)^3*(x(1) - x(2)^4 - 1); 
    J(3) = -1;         
    J(4) = 6*x(1) + 100*(x(1) - x(2)^4 + 1);
    J(5) = -400*x(2)^3*(x(1) - x(2)^4 + 1); 
    J(6) = -1;                              
end

% Lagrangian Hessian (lower triangular)
function H = lagrangian_hessian(x, ~, y)
    H = zeros(3,1);
    H(1) = -106*y(1) - 106*y(2);
    H(2) = 400*x(2)^3*y(1) + 400*x(2)^3*y(2);
    H(3) = 1200*x(2)^2*y(2)*(- x(2)^4 + x(1) + 1) - 1600*x(2)^6*y(2) - 1200*x(2)^2*y(1)*(x(2)^4 - x(1) + 1) - 1600*x(2)^6*y(1);
end

% Lagrangian Hessian operator
function result = lagrangian_hessian_operator(x, ~, y, v)
    h00 = -106*y(1) - 106*y(2);
    h10 = 400*x(2)^3*y(1) + 400*x(2)^3*y(2);
    h11 = 1200*x(2)^2*y(2)*(- x(2)^4 + x(1) + 1) - 1600*x(2)^6*y(2) - 1200*x(2)^2*y(1)*(x(2)^4 - x(1) + 1) - 1600*x(2)^6*y(1);

    result = zeros(3,1);
    result(1) = h00*v(1) + h10*v(2);
    result(2) = h10*v(1) + h11*v(2);
end