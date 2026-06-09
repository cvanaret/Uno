%% UNO high-level interface example
clc, clear

%% HS015 problem
lb = [-inf, -inf]; 
ub = [0.5, inf];
x0 = [-2.0, 1.0];
% options
options = struct();
options.HessianFnc = @hessian_hs015;
% call uno
[x,fval,exitflag,output,lambda,grad,hessian] = uno(@fun_hs015,x0,[],[],[],[],lb,ub,@nlcon_hs015,options);

%% Polak5 problem
x0 = [0.1, 0.1, 0.0];
% options
options = struct();
% options.HessianFnc = @hessian_polak5; % hessian_model='LBFGS' used if HessianFnc is not provided
% call uno
[x,fval,exitflag,output,lambda,grad,hessian] = uno(@fun_polak5,x0,[],[],[],[],[],[],@nlcon_polak5,options);

%% HS015 functions
% Objective function
% objective gradient must be also specified
function [fval,grad] = fun_hs015(x)
    fval = 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
    if nargout > 1
        grad(1) = 400*x(1)^3 - 400*x(1)*x(2) + 2*x(1) - 2;
        grad(2) = 200*(x(2) - x(1)^2);
    end
end

% Constraint function
% constraint gradient must be also specified
function [c,ceq,gradc,gradceq] = nlcon_hs015(x)
    c(1) = -x(1)*x(2)+1;
    c(2) = -(x(1) + x(2)^2);
    ceq = [];
    if nargout > 2
        gradc = [-x(2), -1; ...
                 -x(1), -2*x(2)];
        gradceq = [];
    end
end

% Hessian function (optional)
function H  = hessian_hs015(x,rho,lambda)
    % objective
    H  = [1200*x(1)^2 - 400*x(2) + 2, -400*x(1); 
          -400*x(1), 200] * rho;
    % c(1)
    H  = H + [0, -1;
              -1, 0] * lambda.ineqnonlin(1);
    % c(2)
    H  = H + [0, 0;
              0, 1] * lambda.ineqnonlin(2);
end

%% Polak5 functions
% Objective function
% objective gradient and hessian must be also specified
function [fval,grad] = fun_polak5(x)
    fval = x(3);
    if nargout > 1
        grad = [0; 0; 1];
    end
end

% Constraint function
% constraint gradient and hessian must be also specified
% constraint hessian is specified as a multidimensional matrix
function [c,ceq,gradc,gradceq] = nlcon_polak5(x)
    c(1) = -x(3) + 3*x(1)^2 + 50*(x(1) - x(2)^4 - 1)^2;
    c(2) = -x(3) + 3*x(1)^2 + 50*(x(1) - x(2)^4 + 1)^2;
    ceq = [];
    if nargout > 2
        gradc = [6*x(1) + 100*(x(1) - x(2)^4 - 1), 6*x(1) + 100*(x(1) - x(2)^4 + 1); ...
                 -400*x(2)^3*(x(1) - x(2)^4 - 1), -400*x(2)^3*(x(1) - x(2)^4 + 1); ...
                 -1, -1];
        gradceq = [];
    end
end

% Hessian function (optional)
function H  = hessian_polak5(x,rho,lambda)
    % objective
    H  = zeros(3,3) * rho;
    % c(1)
    H  = H + [106, -400*x(2)^3, 0; 
              -400*x(2)^3, -1200*x(2)^2*(x(1)-x(2)^4-1)+1600*x(2)^6, 0; 
              0, 0, 0] * lambda.ineqnonlin(1);
    % c(2)
    H  = H + [106, -400*x(2)^3, 0; ...
              -400*x(2)^3, -1200*x(2)^2*(x(1)-x(2)^4+1)+1600*x(2)^6, 0;
              0, 0, 0] * lambda.ineqnonlin(2);
end