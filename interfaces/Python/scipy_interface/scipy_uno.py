# Copyright (c) 2026 Francois Gallard
# Licensed under the MIT license. See LICENSE file in the project directory for details.
from numbers import Number
from typing import Iterable, Sequence, Any

import numpy as np
import unopy
from scipy.optimize import LinearConstraint, NonlinearConstraint, OptimizeResult
from scipy.optimize._minimize import standardize_bounds, _validate_bounds, Bounds

Inf = float("inf")
AVAILABLE_METHODS = ["filterslp", "filtersqp", "funnelsqp", "ipopt"]


def minimize(
    fun: callable,
    x0: np.ndarray,
    args: tuple = (),
    method: str = "filtersqp",
    jac: callable = None,
    bounds: Iterable | None | Bounds = None,
    constraints: Iterable[dict[str:Any] | NonlinearConstraint | LinearConstraint] = (),
    tol: float | None = None,
    options: dict | None = None,
) -> OptimizeResult:
    """A scipy.optimize.minimize like interface for Uno.

    Minimization of scalar function of one or more variables under bounds and general constraints.

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
        derivatives ("fun", "jac" functions).
    method : str
        The UNO preset.  Should be one of

        - 'filtersqp'
        - 'filterslp
        - "funnelsqp"
        - "ipopt"

    jac : {callable}, optional
        Method for computing the gradient vector.
        It should be a function that returns the gradient
        vector::

            jac(x, *args) -> array_like, shape (n,)

        where ""x"" is an array with shape (n,) and ""args"" is a tuple with
        the fixed parameters.

    bounds : sequence or "Bounds", optional
        Bounds on variables. There are two ways to specify the bounds:

        1. Instance of "Bounds" class.
        2. Sequence of ""(min, max)"" pairs for each element in "x". None
           is used to specify no bound.

    constraints : {Constraint, dict} or List of {Constraint, dict}, optional
        Constraints definition.

        Available constraints are:

        - "LinearConstraint"
        - "NonlinearConstraint"
        - list of dictionaries.

        For lists of dictionaries, each dictionary with fields:

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
        equal to "tol".

    options : dict, optional
        A dictionary of solver options.

        The available options and their types are:

        "primal_tolerance" : float
        "dual_tolerance" : float
        "loose_primal_tolerance" : float
        "loose_dual_tolerance" : float
        "loose_tolerance_consecutive_iteration_threshold" : int
        "max_iterations" or "max_iter: int
        "time_limit" : float
        "print_solution" : bool
        "unbounded_objective_threshold" : float
        "enforce_linear_constraints" : bool
        "logger" : str,
        "constraint_relaxation_strategy" : str,
        "inequality_handling_method" : str,
        "globalization_mechanism" : str
        "globalization_strategy" : str,
        "hessian_model" : str,
        "inertia_correction_strategy" : str,
        "scale_functions" : bool
        "function_scaling_threshold" : float
        "function_scaling_factor" : float
        "scale_residuals" : bool
        "progress_norm" : str,
        "residual_norm" : str,
        "residual_scaling_threshold" : float
        "protect_actual_reduction_against_roundoff" : bool
        "print_subproblem" : bool
        "armijo_decrease_fraction" : float
        "armijo_tolerance" : float
        "switching_delta" : float
        "switching_infeasibility_exponent" : float
        "filter_type" : str,
        "filter_beta" : float
        "filter_gamma" : float
        "filter_ubd" : float
        "filter_fact" : float
        "filter_capacity" : int
        "filter_sufficient_infeasibility_decrease_factor" : float
        "nonmonotone_filter_number_dominated_entries" : int
        "funnel_kappa" : float
        "funnel_beta" : float
        "funnel_gamma" : float
        "funnel_ubd" : float
        "funnel_fact" : float
        "funnel_update_strategy" : int
        "funnel_require_acceptance_wrt_current_iterate" : bool
        "LS_backtracking_ratio" : float
        "LS_min_step_length" : float
        "LS_scale_duals_with_step_length" : bool
        "regularization_failure_threshold" : float
        "regularization_initial_value" : float
        "regularization_increase_factor" : float
        "primal_regularization_initial_factor" : float
        "dual_regularization_fraction" : float
        "primal_regularization_lb" : float
        "primal_regularization_decrease_factor" : float
        "primal_regularization_fast_increase_factor" : float
        "primal_regularization_slow_increase_factor" : float
        "threshold_unsuccessful_attempts" : int
        "TR_radius" : float
        "TR_increase_factor" : float
        "TR_decrease_factor" : float
        "TR_aggressive_decrease_factor" : float
        "TR_activity_tolerance" : float
        "TR_min_radius" : float
        "TR_radius_reset_threshold" : float
        "switch_to_optimality_requires_linearized_feasibility" : bool
        "l1_constraint_violation_coefficient" : float
        "barrier_initial_parameter" : float
        "barrier_default_multiplier" : float
        "barrier_tau_min" : float
        "barrier_k_sigma" : float
        "barrier_smax" : float
        "barrier_k_mu" : float
        "barrier_theta_mu" : float
        "barrier_k_epsilon" : float
        "barrier_update_fraction" : float
        "barrier_regularization_exponent" : float
        "barrier_small_direction_factor" : float
        "barrier_push_variable_to_interior_k1" : float
        "barrier_push_variable_to_interior_k2" : float
        "barrier_damping_factor" : float
        "least_square_multiplier_max_norm" : float
        "BQPD_kmax_heuristic" : str
        "QP_solver" : str
        "LP_solver" : str
        "linear_solver" : str


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
    Let us consider the problem of minimizing the Rosenbrock function under constraints:
    0.1 <= x**2 <= 0.8.

    >>> res = minimize(rosen, x0 = np.array([1.3, 0.7, 0.8]), jac=rosen_der, method="filtersqp", tol=1e-3,
    >>>                constraints=[NonlinearConstraint(fun=lambda x: x**2,  jac=lambda x: 2 * np.diag(x),
    >>>                                                 lb=np.full(3, 0.1), ub=np.full(3, 0.8))],
    >>>                options={"max_iterations": 10000})
    >>> res.x

    array([ 0.894427, 0.801979, 0.643170])

    To define linear and non-linear constraints, please read the scipy.optimize. minimize documentation:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html

    """
    x0 = np.atleast_1d(np.asarray(x0))

    if x0.ndim != 1:
        raise ValueError("'x0' must only have one dimension.")

    if x0.dtype.kind in np.typecodes["AllInteger"]:
        x0 = np.asarray(x0, dtype=float)

    if not isinstance(args, tuple):
        args = (args,)

    if method not in AVAILABLE_METHODS:
        raise ValueError(f"Unknown solver {method}")

    def _objective(_, x: Sequence, objective_value: Sequence, user_data: Sequence):
        x = np.fromiter(x, dtype="float64")
        if args:
            val = fun(x, *args)
        else:
            val = fun(x)
        if isinstance(val, np.ndarray):
            val = val[0]
        objective_value[0] = val
        return 0

    n_constr = 0
    constr_upper_bounds = []
    constr_lower_bounds = []
    for constr in constraints:
        if isinstance(constr, dict):
            c_type = constr["type"]
            val = constr["fun"](x0)
            if isinstance(val, Number):
                nc = 1
            else:
                nc = len(val)
            if c_type == "ineq":
                ub = np.full(nc, Inf)
                lb = np.zeros(nc)
            elif c_type == "eq":
                ub = np.zeros(nc)
                lb = np.zeros(nc)
            else:
                raise ValueError(f"Unknown constraint type: {c_type}.")
        elif isinstance(constr, LinearConstraint):
            nc = constr.A.shape[0]
            lb = constr.lb
            ub = constr.ub

        elif isinstance(constr, NonlinearConstraint):
            val = constr.fun(x0)
            if isinstance(val, Number):
                nc = 1
            else:
                nc = len(val)

            lb = constr.lb
            ub = constr.ub
        else:
            raise TypeError(f"Unsupported constraint type {type(constr)}")
        n_constr += nc
        if not isinstance(lb, np.ndarray):
            raise ValueError(f"Unsupported lower bound type {type(lb)}")
        if not isinstance(ub, np.ndarray):
            raise ValueError(f"Unsupported upper bound type {type(lb)}")
        constr_lower_bounds.append(lb)
        constr_upper_bounds.append(ub)

    def _constraints(
        _, nc: int, x: Sequence, constraint_values: Sequence, user_data: Sequence
    ):
        x = np.fromiter(x, dtype="float64")
        i_max = 0
        for constr in constraints:
            if isinstance(constr, dict):
                val = constr["fun"](x)
            else:
                if isinstance(constr, LinearConstraint):
                    val = constr.A @ x
                else:
                    val = constr.fun(x)
            if isinstance(val, Number):
                n_v = 1
            else:
                n_v = len(val)
            for i, v in enumerate(val):
                constraint_values[i_max + i] = v
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
            gradient[i] = val[i]
        return 0

    def _constraint_jacobian(_, ncnz, x, jacobian_values, user_data):
        x = np.fromiter(x, dtype="float64")
        i_max = 0
        for constr in constraints:
            if isinstance(constr, dict):
                val = constr["jac"](x)
            else:
                if isinstance(constr, LinearConstraint):
                    val = constr.A
                else:
                    val = constr.jac(x)

            val_flat = val.flatten()
            for i, v in enumerate(val_flat):
                jacobian_values[i_max + i] = v
            i_max += val_flat.size
        return 0

    if bounds is not None:
        # convert to new-style bounds so we only have to consider one case
        _bounds = standardize_bounds(bounds, x0, "new")
        _bounds: Bounds = _validate_bounds(_bounds, x0, method)

        lb = _bounds.lb.tolist()
        ub = _bounds.ub.tolist()
    else:
        lb = [-Inf] * number_variables
        ub = [Inf] * number_variables

    model = unopy.Model(
        unopy.PROBLEM_NONLINEAR, number_variables, lb, ub, unopy.ZERO_BASED_INDEXING
    )
    if jac is not None:
        model.set_objective(unopy.MINIMIZE, _objective, _objective_gradient)
    else:
        raise ValueError("Jacobian of objective is required.")

    if n_constr > 0:
        number_jacobian_nonzeros = number_variables * n_constr
        jacobian_row_indices = np.tile(np.arange(number_variables), n_constr).tolist()
        jacobian_column_indices = np.repeat(
            np.arange(number_variables), n_constr
        ).tolist()
        constr_lower_bounds = np.concatenate(constr_lower_bounds)
        constr_upper_bounds = np.concatenate(constr_upper_bounds)

        model.set_constraints(
            n_constr,
            _constraints,
            constr_lower_bounds,
            constr_upper_bounds,
            number_jacobian_nonzeros,
            jacobian_row_indices,
            jacobian_column_indices,
            _constraint_jacobian,
        )
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

    if method in ["filterslp"]:
        if options is None or "LP_solver" not in options:
            uno_solver.set_option("LP_solver", "HiGHS")
    if method in ["ipopt"]:
        if options is None or "linear_solver" not in options:
            try:
                uno_solver.set_option("linear_solver", "MUMPS")
            except ValueError:
                uno_solver.set_option("linear_solver", "MA57")

    if options is not None:
        if "max_iter" in options:
            uno_solver.set_option("max_iter", options["max_iter"])
        for opt_key, opt_value in options.items():
            uno_solver.set_option(opt_key, opt_value)

    result = uno_solver.optimize(model)
    return OptimizeResult(
        x=result.primal_solution,
        success=str(result.optimization_status) == "OptimizationStatus.SUCCESS",
        status=result.solution_status,
        message=result.optimization_status,
        fun=result.solution_objective,
        nfev=result.number_objective_evaluations,
        njev=result.number_jacobian_evaluations,
        nhev=result.number_hessian_evaluations,
        nit=result.number_iterations,
        maxcv=result.solution_primal_feasibility,
    )
