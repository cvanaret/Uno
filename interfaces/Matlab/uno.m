function [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = uno(FUN,X,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
% UNO - UNO high-level interface
%   Solve an optimization model described using UNO, a modern and modular solver for nonlinearly constrained optimization.
%
%   Syntax
%     x = uno(fun,x0)
%     x = uno(fun,x0,A,b)
%     x = uno(fun,x0,A,b,Aeq,beq)
%     x = uno(fun,x0,A,b,Aeq,beq,lb,ub)
%     x = uno(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
%     x = uno(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
%     options = uno('defaults')
%     options = uno('defaults', preset)
%
%   Input Arguments
%     fun - Function to minimize
%       function handle | function name
%     x0 - Initial point
%       real vector | real array
%     A - Linear inequality constraints
%       real matrix
%     b - Linear inequality constraints
%       real vector
%     Aeq - Linear equality constraints
%       real matrix
%     beq - Linear equality constraints
%       real vector
%     lb - Lower bounds
%       real vector | real array
%     ub - Upper bounds
%       real vector | real array
%     nonlcon - Nonlinear constraints
%       function handle | function name
%     options - UNO options, as returned by <a href="matlab:eval('uno(''defaults'')')">uno('defaults')</a>  
%       structure
%     preset - Option preset
%       character vector | string
%
%   Output Arguments
%     x - Solution
%       real vector | real array
%     fval - Objective function value at solution
%       real number
%     exitflag - Reason fmincon stopped
%       integer
%     output - Information about the optimization process
%       structure
%     lambda - Lagrange multipliers at the solution
%       structure
%     grad - Gradient at the solution
%       real vector
%     hessian - Hessian matrix
%       real matrix
%     options - UNO options
%       structure
%
%   Examples
%     <a href="matlab:open([fileparts(which('uno')) filesep 'example/example_uno.m'])">example_uno</a>
%
% Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
% Licensed under the MIT license. See LICENSE file in the project directory for details.

% return the options if 'defaults' and preset
if nargin<=2 && nargout <= 1 && strcmpi(FUN,'defaults')
    if nargin==1
        X = uno_options();
    else
        X = uno_options(X);
    end
    return
end

% default args
if nargin < 10
    options = struct();
    if nargin < 9
        NONLCON = [];
        if nargin < 8
            UB = [];
            if nargin < 7
                LB = [];
                if nargin < 6
                    Beq = [];
                    if nargin < 5
                        Aeq = [];
                        if nargin < 4
                            B = [];
                            if nargin < 3
                                A = [];
                            end
                        end
                    end
                end
            end
        end
    end
end

% Column vector
X = X(:);
Beq = Beq(:);
B = B(:);

% Check for complex X0
if ~isreal(X)
    error('uno:ComplexX0', ...
        getString(message('optimlib:commonMsgs:ComplexX0','Uno')));
end

% Check for empty X
if isempty(X)
    error(message('optimlib:fmincon:EmptyX'));
end

% Set empty linear constraints
if isempty(Aeq)
    Aeq = zeros(0,length(X));
end
if isempty(A)
    A = zeros(0,length(X));
end
if isempty(Beq)
    Beq = zeros(0,1);
end
if isempty(B)
    B = zeros(0,1);
end

% Check linear constraint sizes
if size(Aeq,2) ~= length(X)
    error(message('optimlib:fmincon:WrongNumberOfColumnsInAeq', length(X)))
end
if size(Aeq,1) ~= length(Beq)
    error(message('optimlib:fmincon:AeqAndBeqInconsistent'))
end
if size(A,2) ~= length(X)
    error(message('optimlib:fmincon:WrongNumberOfColumnsInA', length(X)))
end
if size(A,1) ~= length(B)
    error(message('optimlib:fmincon:AeqAndBinInconsistent'))
end

% Check the bounds
[X,l,u,msg] = checkbounds(X,LB,UB,length(X));
if ~isempty(msg)
    EXITFLAG = -2;
    [FVAL,LAMBDA,GRAD,HESSIAN] = deal([]);
    OUTPUT.iterations = 0;
    OUTPUT.funcCount = 0;
    OUTPUT.constrviolation  = [];
    OUTPUT.firstorderopt = [];
    OUTPUT.message = msg;
    return
end

% Check [fval,grad,hess] = FUN(x)
if ~isa(FUN,'function_handle')
    error('FUN must be a function handle.');
end
if (nargin(FUN)~=1) || (nargout(FUN)<3)
    error('FUN must have 1 input and 3 output argument(s).')
end
try
    [fval0,grad0,hess0] = FUN(X);
catch
    error('Failure in objective function evaluation.')
end
if ~isnumeric(fval0) || ~isnumeric(grad0) || ~isnumeric(hess0)
    error('FUM must return numeric values.')
end
if ~isscalar(fval0)
    error('FUM must return a scalar value.')
end
if length(grad0) ~= length(X)
    error('Dimension of grad is inconsistent with the length of x0.')
end
if ~isempty(hess0) && ~isequal(size(hess0), [length(X), length(X)])
    error('Dimension of hess is inconsistent with the length of x0.')
end

% Check [c,ceq,gradc,gradceq,hessc,hessceq] = NONLCON(x)
if isempty(NONLCON)
    NONLCON = @(x) deal([]);
else
    if ~isa(NONLCON,'function_handle')
        error('NLCON must be a function handle.')
    end
    if (nargin(NONLCON)~=1) || (nargout(NONLCON)<6)
        error('NONLCON must have 1 input and 6 output argument(s).')
    end
end
try
    [c0,ceq0,gradc0,gradceq0,hessc0,hessceq0] = NONLCON(X);
catch
    error('Failure in nonlinear constraint function evaluation.')
end
if ~isnumeric(c0) || ~isnumeric(ceq0) || ~isnumeric(gradc0) || ~isnumeric(gradceq0) ...
         || ~isnumeric(hessc0) || ~isnumeric(hessceq0)
    error('NONLCON must return numeric values.')
end
if length(c0) ~= size(gradc0,2)
    error('Column dimension of gradc is inconsistent with the length of c.')
end
if length(ceq0) ~= size(gradceq0,2)
    error('Column dimension of gradceq is inconsistent with the length of ceq.')
end
if ~isempty(c0) && length(X) ~= size(gradc0,1)
    error('Row dimension of gradc is inconsistent with the length of x0.')
end
if ~isempty(ceq0) && length(X) ~= size(gradceq0,1)
    error('Row dimension of gradceq is inconsistent with the length of x0.')
end
if ~isempty(c0) && ~isempty(hessc0) && ~isequal(size(hessc0), [length(X), length(X), length(c0)])
    error('Dimension of hess is inconsistent with the length of x0 and c.')
end
if ~isempty(ceq0) && ~isempty(hessceq0) && ~isequal(size(hessceq0), [length(X), length(X), length(ceq0)])
    error('Dimension of hess is inconsistent with the length of x0 and ceq.')
end

% Constraint indexes
ind_eqlin = 1 : numel(Beq);
ind_ineqlin = (numel(Beq)+1) : (numel(Beq)+numel(B));
ind_ineqnonlin = (numel(Beq)+numel(B)+1) : (numel(Beq)+numel(B)+numel(c0));
ind_eqnonlin = (numel(Beq)+numel(B)+numel(c0)+1) : (numel(Beq)+numel(B)+numel(c0)+numel(ceq0));

% Create the optimization problem
model.problem_type = 'N';
model.base_indexing = 1;
model.number_variables = length(X);
model.variables_lower_bounds = l;
model.variables_upper_bounds = u;
model.optimization_sense = 1;
model.objective_function = @(x) objective_function(x, FUN); 
model.objective_gradient = @(x) objective_gradient(x, FUN);
model.number_constraints = length(Beq) + length(B) + length(c0) + length(ceq0);
model.constraints_lower_bounds = [zeros(size(Beq)); -inf(size(B)); -inf(length(c0),1); zeros(length(ceq0),1)];
model.constraints_upper_bounds = zeros(model.number_constraints,1);
model.constraint_function = @(x) constraint_function(x,Aeq,Beq,A,B,NONLCON);
[irj, jcj] = find(ones(model.number_constraints,length(X)));
model.jacobian_row_indices = irj;
model.jacobian_column_indices = jcj;
model.number_jacobian_nonzeros = length(irj);
model.constraint_jacobian = @(x) constraint_jacobian(x,Aeq,Beq,A,B,NONLCON,irj,jcj);
model.lagrangian_sign_convention = 1;
model.hessian_triangular_part = 'L';
[irh, jch] = find(tril(ones(length(X),length(X))));
model.hessian_row_indices = irh(:);
model.hessian_column_indices = jch(:);
model.number_hessian_nonzeros =  length(irh(:));
model.lagrangian_hessian = @(x,rho,y) lagrangian_hessian(x,rho,y,FUN,NONLCON,ind_ineqnonlin,ind_eqnonlin,irh,jch);
model.initial_primal_iterate = X;

% Call UNO
result = uno_optimize(model, options);

% Results
X = result.primal_solution(1:length(X));
FVAL = result.solution_objective;
switch result.optimization_status
    case 0
        switch result.solution_status
            case {0,3,5,6}
                EXITFLAG = -2;
            case {1,2}
                EXITFLAG = 1;
            case 4
                EXITFLAG = 2;
        end
    case 1
        EXITFLAG = 0;
    case {2,3,4}
        EXITFLAG = -2;
    case 5
        EXITFLAG = -1;
end
OUTPUT.iterations = result.number_iterations;
OUTPUT.funcCount = result.number_objective_evaluations;
OUTPUT.constrviolation  = result.solution_primal_feasibility;
OUTPUT.firstorderopt = result.solution_stationarity;
LAMBDA.lower = result.lower_bound_dual_solution(1:length(X));
LAMBDA.upper = result.upper_bound_dual_solution(1:length(X));
LAMBDA.eqlin = result.constraint_dual_solution(ind_eqlin);
LAMBDA.ineqlin = result.constraint_dual_solution(ind_ineqlin);
LAMBDA.ineqnonlin = result.constraint_dual_solution(ind_ineqnonlin);
LAMBDA.eqnonlin = result.constraint_dual_solution(ind_eqnonlin);
[~,GRAD] = FUN(X);
HESSIAN = hessian(X, 1, LAMBDA.ineqnonlin, LAMBDA.eqnonlin, FUN, NONLCON);

end

% Model functions
function fval = objective_function(x,FUN)
    fval = FUN(x);
end

function grad = objective_gradient(x,FUN)
    [~, grad] = FUN(x);
end

function constr = constraint_function(x,Aeq,Beq,A,B,NONLCON)
    [c,ceq] = NONLCON(x);
    constr = [Aeq*x-Beq; ...
              A*x-B; ...
              c(:); ...
              ceq(:)];   
end

function J = constraint_jacobian(x,Aeq,~,A,~,NONLCON,ir,jc)
    [~,~,gradc,gradceq] = NONLCON(x);
    J = [Aeq; ...
         A; ...
         gradc'; ...
         gradceq'];
    J = J(ir + (jc-1)*size(J,1));
end

function H = hessian(x, rho, yneq, yeq, FUN, NONLCON)
    [~, ~, hess] = FUN(x);
    [~,~,~,~,hessc,hessceq] = NONLCON(x);
    H = zeros(length(x), length(x));
    if ~isempty(hess)
        H = H + rho*hess;
    end
    if ~isempty(hessc)
        H = H + sum(hessc .* reshape(yneq, [1,1,length(yneq)]), 3);
    end
    if ~isempty(hessceq)
        H = H + sum(hessceq .* reshape(yeq, [1,1,length(yeq)]), 3);
    end
end
 
function H = lagrangian_hessian(x, rho, y, FUN, NONLCON, icneq, iceq, ir, jc)
    H = hessian(x, rho, y(icneq), y(iceq), FUN, NONLCON);
    H = H(ir + (jc-1)*size(H,1));
end