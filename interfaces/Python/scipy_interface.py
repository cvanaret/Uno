from typing import Iterable

import numpy as np
import unopy
from scipy.optimize import LinearConstraint, NonlinearConstraint, OptimizeResult
from scipy.optimize._minimize import standardize_bounds, _validate_bounds, Bounds

Inf = float("inf")


def minimize(fun: callable, x0: np.ndarray, args: tuple = (), method: str = "filtersqp", jac: callable = None,
             hess: callable = None,
             hessp: callable = None, bounds: Iterable | None | Bounds = None,
             constraints: Iterable = (),
             tol: float | None = None, callback=None, options: dict | None = None) -> OptimizeResult:
    """Minimization of scalar function of one or more variables.

    Parameters
    ----------
    fun : callable
        The objective function to be minimized::

            fun(x, *args) -> float

        where ""x"" is a 1-D array with shape (n,) and ""args""
        is a tuple of the fixed parameters needed to completely
        specify the function.

        Suppose the callable has signature ""f0(x, *my_args, **my_kwargs)"", where
        ""my_args"" and ""my_kwargs"" are required positional and keyword arguments.
        Rather than passing ""f0"" as the callable, wrap it to accept
        only ""x""; e.g., pass ""fun=lambda x: f0(x, *my_args, **my_kwargs)"" as the
        callable, where ""my_args"" (tuple) and ""my_kwargs"" (dict) have been
        gathered before invoking this function.
    x0 : ndarray, shape (n,)
        Initial guess. Array of real elements of size (n,),
        where ""n"" is the number of independent variables.
    args : tuple, optional
        Extra arguments passed to the objective function and its
        derivatives ("fun", "jac" and "hess" functions).
    method : str or callable, optional
        Type of solver.  Should be one of

        - 'filtersqp'
        - 'filterslp'

    jac : {callable}, optional
        Method for computing the gradient vector.
        It should be a function that returns the gradient
        vector::

            jac(x, *args) -> array_like, shape (n,)

        where ""x"" is an array with shape (n,) and ""args"" is a tuple with
        the fixed parameters. If "jac" is a Boolean and is True, "fun" is
        assumed to return a tuple ""(f, g)"" containing the objective
        function and the gradient.
        Methods 'Newton-CG', 'trust-ncg', 'dogleg', 'trust-exact', and
        'trust-krylov' require that either a callable be supplied, or that
        "fun" return the objective and gradient.
        If None or False, the gradient will be estimated using 2-point finite
        difference estimation with an absolute step size.
        Alternatively, the keywords  {'2-point', '3-point', 'cs'} can be used
        to select a finite difference scheme for numerical estimation of the
        gradient with a relative step size. These finite difference schemes
        obey any specified "bounds".

    hessp : callable, optional
        Hessian of objective function times an arbitrary vector p. Only for
        Newton-CG, trust-ncg, trust-krylov, trust-constr.
        Only one of "hessp" or "hess" needs to be given. If "hess" is
        provided, then "hessp" will be ignored. "hessp" must compute the
        Hessian times an arbitrary vector::

            hessp(x, p, *args) ->  ndarray shape (n,)

        where ""x"" is a (n,) ndarray, ""p"" is an arbitrary vector with
        dimension (n,) and ""args"" is a tuple with the fixed
        parameters.
    bounds : sequence or "Bounds", optional
        Bounds on variables for Nelder-Mead, L-BFGS-B, TNC, SLSQP, Powell,
        trust-constr, COBYLA, and COBYQA methods. There are two ways to specify
        the bounds:

        1. Instance of "Bounds" class.
        2. Sequence of ""(min, max)"" pairs for each element in "x". None
           is used to specify no bound.

    constraints : {Constraint, dict} or List of {Constraint, dict}, optional
        Constraints definition.

        Constraints for 'trust-constr', 'cobyqa', and 'cobyla' are defined as a single
        object or a list of objects specifying constraints to the optimization problem.
        Available constraints are:

        - "LinearConstraint"
        - "NonlinearConstraint"

        Constraints for COBYLA, SLSQP are defined as a list of dictionaries.
        Each dictionary with fields:

        type : str
            Constraint type: 'eq' for equality, 'ineq' for inequality.
        fun : callable
            The function defining the constraint.
        jac : callable, optional
            The Jacobian of "fun" (only for SLSQP).
        args : sequence, optional
            Extra arguments to be passed to the function and Jacobian.

        Equality constraint means that the constraint function result is to
        be zero whereas inequality means that it is to be non-negative.

    tol : float, optional
        Tolerance for termination. When "tol" is specified, the selected
        minimization algorithm sets some relevant solver-specific tolerance(s)
        equal to "tol". For detailed control, use solver-specific
        options.

    options : dict, optional
        A dictionary of solver options.

        maxiter : int
            Maximum number of iterations to perform. Depending on the
            method each iteration may use several function evaluations.

            For "TNC" use "maxfun" instead of "maxiter".
        disp : bool
            Set to True to print convergence messages.

        For method-specific options, see :func:"show_options()".

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ""OptimizeResult"" object.
        Important attributes are: ""x"" the solution array, ""success"" a
        Boolean flag indicating if the optimizer exited successfully and
        ""message"" which describes the cause of the termination. See
        "OptimizeResult" for a description of other attributes.

    Examples
    --------
    Let us consider the problem of minimizing the Rosenbrock function. This
    function (and its respective derivatives) is implemented in "rosen"
    (resp. "rosen_der", "rosen_hess") in the "scipy.optimize".

    >>> from unopy import minimize

    A simple application of the *Nelder-Mead* method is:

    >>> x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
    >>> res = minimize(rosen, x0, method='filtersqp', tol=1e-6)
    >>> res.x
    array([ 1.,  1.,  1.,  1.,  1.])

    Now using the *BFGS* algorithm, using the first derivative and a few
    options:

    >>> res = minimize(rosen, x0, method='BFGS', jac=rosen_der,
    ...                options={'gtol': 1e-6, 'disp': True})
    Optimization terminated successfully.
             Current function value: 0.000000
             Iterations: 26
             Function evaluations: 31
             Gradient evaluations: 31
    >>> res.x
    array([ 1.,  1.,  1.,  1.,  1.])
    >>> print(res.message)
    Optimization terminated successfully.

    Next, consider a minimization problem with several constraints (namely
    Example 16.4 from [5]_). The objective function is:

    >>> fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2

    There are three constraints defined as:

    >>> cons = ({'type': 'ineq', 'fun': lambda x:  x[0] - 2 * x[1] + 2},
    ...         {'type': 'ineq', 'fun': lambda x: -x[0] - 2 * x[1] + 6},
    ...         {'type': 'ineq', 'fun': lambda x: -x[0] + 2 * x[1] + 2})

    And variables must be positive, hence the following bounds:

    >>> bnds = ((0, None), (0, None))

    The optimization problem is solved using the SLSQP method as:

    >>> res = minimize(fun, (2, 0), method='SLSQP', bounds=bnds, constraints=cons)

    It should converge to the theoretical solution ""[1.4 ,1.7]"". *SLSQP* also
    returns the multipliers that are used in the solution of the problem. These
    multipliers, when the problem constraints are linear, can be thought of as the
    Karush-Kuhn-Tucker (KKT) multipliers, which are a generalization
    of Lagrange multipliers to inequality-constrained optimization problems ([20]_).

    Notice that at the solution, the first constraint is active. Let's evaluate the
    function at solution:

    >>> cons[0]['fun'](res.x)
    np.float64(1.4901224698604665e-09)

    Also, notice that at optimality there is a non-zero multiplier:

    >>> res.multipliers
    array([0.8, 0. , 0. ])

    This can be understood as the local sensitivity of the optimal value of the
    objective function with respect to changes in the first constraint. If we
    tighten the constraint by a small amount ""eps"":

    >>> eps = 0.01
    >>> cons[0]['fun'] = lambda x: x[0] - 2 * x[1] + 2 - eps

    we expect the optimal value of the objective function to increase by
    approximately ""eps * res.multipliers[0]"":

    >>> eps * res.multipliers[0]  # Expected change in f0
    np.float64(0.008000000027153205)
    >>> f0 = res.fun  # Keep track of the previous optimal value
    >>> res = minimize(fun, (2, 0), method='SLSQP', bounds=bnds, constraints=cons)
    >>> f1 = res.fun  # New optimal value
    >>> f1 - f0
    np.float64(0.008019998807885509)

    """
    x0 = np.atleast_1d(np.asarray(x0))

    if x0.ndim != 1:
        raise ValueError("'x0' must only have one dimension.")

    if x0.dtype.kind in np.typecodes["AllInteger"]:
        x0 = np.asarray(x0, dtype=float)

    if not isinstance(args, tuple):
        args = (args,)

    if method not in ["filterslp", "filtersqp", "funnelsqp", "ipopt"]:
        raise ValueError(f'Unknown solver {method}')

    def _objective(_, x, objective_value, user_data):
        x= np.fromiter(x, dtype="float64")
        if args:
            val = fun(x, *args)
        else:
            val = fun(x)
        if isinstance(val, np.ndarray):
            val=val[0]
        objective_value[0]=val
        return 0

    n_constr = 0
    constr_upper_bounds = []
    constr_lower_bounds = []
    for constr in constraints:
        if isinstance(constr, dict):
            c_type = constr["type"]
            val = constr["fun"](x0).size
            if isinstance(val, float):
                n_constr += 1
            else:
                n_constr += len(val)
            if c_type == "ineq":
                constr_upper_bounds.append(Inf)
                constr_lower_bounds.append(0.)
            elif c_type == "eq":
                constr_upper_bounds.append(0.)
                constr_lower_bounds.append(0.)
            else:
                raise ValueError(f"Unknown constraint type: {c_type}.")
        elif isinstance(constr, LinearConstraint):
            constr_lower_bounds.append(constr.lb)
            constr_upper_bounds.append(constr.ub)
            n_constr += constr.A.shape[0]
        elif isinstance(constr, NonlinearConstraint):
            val = constr.fun(x0).size
            if isinstance(val, float):
                n_constr += 1
            else:
                n_constr += len(val)

            constr_lower_bounds.append(constr.lb)
            constr_upper_bounds.append(constr.ub)

    def _constraints(_, nc, x, constraint_values, user_data):
        x = np.fromiter(x, dtype="float64")
        i_max = 0
        for constr in constraints:
            if isinstance(constr, dict):
                val = constr['fun'](x)
            else:
                val=constr.fun(x)
            if isinstance(val, float):
                n_v = 1
            else:
                n_v = len(val)
            constraint_values[i_max:i_max + n_v] = val
            i_max += n_v
        return 0

    number_variables = x0.size
    def _objective_gradient(_, x, gradient, user_data):
        x = np.fromiter(x, dtype="float64")
        if args:
            val = jac(x, *args)
        else:
            val = jac(x)
        for i in range(number_variables):
            gradient[i]=val[i]
        return 0



    def _constraint_jacobian(_, ncnz, x, jacobian_values, user_data):
        x = np.fromiter(x, dtype="float64")
        i_max = 0
        for constr in constraints:
            if isinstance(constr, dict):
                val = constr["jac"](x)
            else:
                val = constr.jac(x)
            if isinstance(val, float):
                n_v = number_variables
            else:
                n_v = len(val)*number_variables
            jacobian_values[i_max:i_max + n_v] = val.flatten()
            i_max += n_v
        return 0


    if bounds is not None:
        # convert to new-style bounds so we only have to consider one case
        _bounds = standardize_bounds(bounds, x0, 'new')
        _bounds: Bounds = _validate_bounds(_bounds, x0, method)

        lb = _bounds.lb.tolist()
        ub = _bounds.ub.tolist()
    else:
        lb = [-Inf] * number_variables
        ub = [Inf] * number_variables

    model = unopy.Model(unopy.PROBLEM_NONLINEAR, number_variables, lb, ub,
                        unopy.ZERO_BASED_INDEXING)
    if jac is not None:
        model.set_objective(unopy.MINIMIZE, _objective, _objective_gradient)
    else:
        raise ValueError("Jacobian of objective is required.")

    if n_constr>0:
        number_jacobian_nonzeros = number_variables * n_constr
        jacobian_row_indices = np.tile(np.arange(number_variables), n_constr).tolist()  # [0, 1, 0, 1]
        print("jacobian_row_indices", jacobian_row_indices)
        jacobian_column_indices = np.repeat(np.arange(number_variables), n_constr).tolist()  # [0, 0, 1, 1]
        print("jacobian_column_indices", jacobian_column_indices)
        model.set_constraints(n_constr, _constraints, constr_lower_bounds, constr_upper_bounds,
                              number_jacobian_nonzeros, jacobian_row_indices, jacobian_column_indices,
                              _constraint_jacobian)
    # model.set_lagrangian_hessian(number_hessian_nonzeros, hessian_triangular_part, hessian_row_indices,
    #                             hessian_column_indices, _lagrangian_hessian, unopy.MULTIPLIER_NEGATIVE)
    # model.set_lagrangian_hessian_operator(lagrangian_hessian_operator, unopy.MULTIPLIER_NEGATIVE)
    model.set_initial_primal_iterate(x0.tolist())
    if args:
        model.set_user_data(args)

    uno_solver = unopy.UnoSolver()

    uno_solver.set_preset(method)
    if tol is not None:
        uno_solver.set_option("primal_tolerance", tol)
        uno_solver.set_option("dual_tolerance", tol)
        uno_solver.set_option("loose_dual_tolerance", tol)
        uno_solver.set_option("loose_tolerance_consecutive_iteration_threshold", tol)
        uno_solver.set_option("primal_tolerance", tol)

    if method in ["filterslp" ]:
        if options is None or "LP_solver" not in options:
            uno_solver.set_option("LP_solver", "HiGHS")
    if method in ["ipopt"]:
        if options is None or "linear_solver" not in options:
            try:
                uno_solver.set_option("linear_solver", "MUMPS")
            except ValueError:
                uno_solver.set_option("linear_solver", "MA57")

    if options is not None :
        if "max_iter" in options:
            uno_solver.set_option("max_iter", options["max_iter"])
        for opt_key, opt_value in options.items():
            uno_solver.set_option(opt_key, opt_value)


    result = uno_solver.optimize(model)
    return OptimizeResult(x=result.primal_solution, success=result.optimization_status == "UNO_SUCCESS",
                   status=result.solution_status, message=result.optimization_status, fun=result.solution_objective,
                   nfev=result.number_objective_evaluations, njev=result.number_jacobian_evaluations,
                   nhev=result.number_hessian_evaluations, nit=result.number_iterations,
                   maxcv=result.solution_primal_feasibility)

#
# "primal_tolerance", OptionType::DOUBLE},
# "dual_tolerance", OptionType::DOUBLE},
# "loose_primal_tolerance", OptionType::DOUBLE},
# "loose_dual_tolerance", OptionType::DOUBLE},
# "loose_tolerance_consecutive_iteration_threshold", OptionType::INTEGER},
# "max_iterations", OptionType::INTEGER},
# "time_limit", OptionType::DOUBLE},
# "print_solution", OptionType::BOOL},
# "unbounded_objective_threshold", OptionType::DOUBLE},
# "enforce_linear_constraints", OptionType::BOOL},
# "logger", OptionType::STRING},
# "constraint_relaxation_strategy", OptionType::STRING},
# "inequality_handling_method", OptionType::STRING},
# "globalization_mechanism",OptionType::STRING},
# "globalization_strategy", OptionType::STRING},
# "hessian_model", OptionType::STRING},
# "inertia_correction_strategy", OptionType::STRING},
# "scale_functions", OptionType::BOOL},
# "function_scaling_threshold", OptionType::DOUBLE},
# "function_scaling_factor", OptionType::DOUBLE},
# "scale_residuals", OptionType::BOOL},
# "progress_norm", OptionType::STRING},
# "residual_norm", OptionType::STRING},
# "residual_scaling_threshold", OptionType::DOUBLE},
# "protect_actual_reduction_against_roundoff", OptionType::BOOL},
# "print_subproblem", OptionType::BOOL},
# "armijo_decrease_fraction", OptionType::DOUBLE},
# "armijo_tolerance", OptionType::DOUBLE},
# "switching_delta", OptionType::DOUBLE},
# "switching_infeasibility_exponent", OptionType::DOUBLE},
# "filter_type", OptionType::STRING},
# "filter_beta", OptionType::DOUBLE},
# "filter_gamma", OptionType::DOUBLE},
# "filter_ubd", OptionType::DOUBLE},
# "filter_fact", OptionType::DOUBLE},
# "filter_capacity", OptionType::INTEGER},
# "filter_sufficient_infeasibility_decrease_factor", OptionType::DOUBLE},
# "nonmonotone_filter_number_dominated_entries", OptionType::INTEGER},
# "funnel_kappa", OptionType::DOUBLE},
# "funnel_beta", OptionType::DOUBLE},
# "funnel_gamma", OptionType::DOUBLE},
# "funnel_ubd", OptionType::DOUBLE},
# "funnel_fact", OptionType::DOUBLE},
# "funnel_update_strategy", OptionType::INTEGER},
# "funnel_require_acceptance_wrt_current_iterate", OptionType::BOOL},
# "LS_backtracking_ratio", OptionType::DOUBLE},
# "LS_min_step_length", OptionType::DOUBLE},
# "LS_scale_duals_with_step_length", OptionType::BOOL},
# "regularization_failure_threshold", OptionType::DOUBLE},
# "regularization_initial_value", OptionType::DOUBLE},
# "regularization_increase_factor", OptionType::DOUBLE},
# "primal_regularization_initial_factor", OptionType::DOUBLE},
# "dual_regularization_fraction", OptionType::DOUBLE},
# "primal_regularization_lb", OptionType::DOUBLE},
# "primal_regularization_decrease_factor", OptionType::DOUBLE},
# "primal_regularization_fast_increase_factor", OptionType::DOUBLE},
# "primal_regularization_slow_increase_factor", OptionType::DOUBLE},
# "threshold_unsuccessful_attempts", OptionType::INTEGER},
# "TR_radius", OptionType::DOUBLE},
# "TR_increase_factor", OptionType::DOUBLE},
# "TR_decrease_factor", OptionType::DOUBLE},
# "TR_aggressive_decrease_factor", OptionType::DOUBLE},
# "TR_activity_tolerance", OptionType::DOUBLE},
# "TR_min_radius", OptionType::DOUBLE},
# "TR_radius_reset_threshold", OptionType::DOUBLE},
# "switch_to_optimality_requires_linearized_feasibility", OptionType::BOOL},
# "l1_constraint_violation_coefficient", OptionType::DOUBLE},
# "barrier_initial_parameter", OptionType::DOUBLE},
# "barrier_default_multiplier", OptionType::DOUBLE},
# "barrier_tau_min", OptionType::DOUBLE},
# "barrier_k_sigma", OptionType::DOUBLE},
# "barrier_smax", OptionType::DOUBLE},
# "barrier_k_mu", OptionType::DOUBLE},
# "barrier_theta_mu", OptionType::DOUBLE},
# "barrier_k_epsilon", OptionType::DOUBLE},
# "barrier_update_fraction", OptionType::DOUBLE},
# "barrier_regularization_exponent", OptionType::DOUBLE},
# "barrier_small_direction_factor", OptionType::DOUBLE},
# "barrier_push_variable_to_interior_k1", OptionType::DOUBLE},
# "barrier_push_variable_to_interior_k2", OptionType::DOUBLE},
# "barrier_damping_factor",
# "least_square_multiplier_max_norm",
# "BQPD_kmax_heuristic",
# "QP_solver"
# "LP_solver"
# "linear_solver"
