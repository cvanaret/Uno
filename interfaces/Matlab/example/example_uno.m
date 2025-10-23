%% UNO high-level interface example
clc, clear

%% HS015 problem
lb = [-inf, -inf]; 
ub = [0.5, inf];
x0 = [-2.0, 1.0];

[x,fval,exitflag,output,lambda,grad,hessian] = uno(@fun_hs015,x0,[],[],[],[],lb,ub,@nlcon_hs015);

%% Polak5 problem
x0 = [0.1, 0.1, 0.0];

[x,fval,exitflag,output,lambda,grad,hessian] = uno(@fun_polak5,x0,[],[],[],[],[],[],@nlcon_polak5);

%% HS015 functions
% Objective function
% objective gradient and hessian must be also specified
function [fval,grad,hess] = fun_hs015(x)
    fval = 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
    if nargout > 1
        grad(1) = 400*x(1)^3 - 400*x(1)*x(2) + 2*x(1) - 2;
        grad(2) = 200*(x(2) - x(1)^2);
    end
    if nargout > 2
        hess = [1200*x(1)^2 - 400*x(2) + 2, -400*x(1); ...
                -400*x(1), 200];
    end
end

% Constraint function
% constraint gradient and hessian must be also specified
% constraint hessian is specified as a multidimensional matrix
function [c,ceq,gradc,gradceq,hessc,hessceq] = nlcon_hs015(x)
    c(1) = -x(1)*x(2)+1;
    c(2) = -(x(1) + x(2)^2);
    ceq = [];
    if nargout > 2
        gradc = [-x(2), -1; ...
                 -x(1), -2*x(2)];
        gradceq = [];
    end
    if nargout > 4
        hessc(:,:,1) = [0, -1; ...
                        -1, 0];
        hessc(:,:,2) = [0, 0; ...
                        0, 1];
        hessceq = [];
    end
end

%% Polak5 functions
% Objective function
% objective gradient and hessian must be also specified
function [fval,grad,hess] = fun_polak5(x)
    fval = x(3);
    if nargout > 1
        grad = [0; 0; 1];
    end
    if nargout > 2
        hess = zeros(3,3);
    end
end

% Constraint function
% constraint gradient and hessian must be also specified
% constraint hessian is specified as a multidimensional matrix
function [c,ceq,gradc,gradceq,hessc,hessceq] = nlcon_polak5(x)
    c(1) = -x(3) + 3*x(1)^2 + 50*(x(1) - x(2)^4 - 1)^2;
    c(2) = -x(3) + 3*x(1)^2 + 50*(x(1) - x(2)^4 + 1)^2;
    ceq = [];
    if nargout > 2
        gradc = [6*x(1) + 100*(x(1) - x(2)^4 - 1), 6*x(1) + 100*(x(1) - x(2)^4 + 1); ...
                 -400*x(2)^3*(x(1) - x(2)^4 - 1), -400*x(2)^3*(x(1) - x(2)^4 + 1); ...
                 -1, -1];
        gradceq = [];
    end
    if nargout > 4
        hessc(:,:,1) = [106, -400*x(2)^3, 0; ...
                       -400*x(2)^3, -1200*x(2)^2*(x(1)-x(2)^4-1)+1600*x(2)^6, 0; ...
                       0, 0, 0];
        hessc(:,:,2) = [106, -400*x(2)^3, 0; ...
                       -400*x(2)^3, -1200*x(2)^2*(x(1)-x(2)^4+1)+1600*x(2)^6, 0;
                        0, 0, 0];
        hessceq = [];
    end
end
