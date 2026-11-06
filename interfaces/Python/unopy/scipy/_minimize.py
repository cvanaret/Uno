# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""A ``scipy.optimize.minimize``-like interface to Uno."""

from __future__ import annotations

from typing import Any
from typing import Callable
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import BFGS
from scipy.optimize import Bounds
from scipy.optimize import HessianUpdateStrategy
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint
from scipy.optimize import OptimizeResult
from scipy.optimize import SR1

import unopy
from unopy.scipy._bounds import _get_variable_bounds
from unopy.scipy._callback import STOP_ITERATION_MESSAGE
from unopy.scipy._callback import _IterationCallback
from unopy.scipy._constraints import _MergedConstraints
from unopy.scipy._constraints import _normalize_constraints
from unopy.scipy._finite_differences import _is_finite_difference_method
from unopy.scipy._guard import _EvaluationErrorGuard
from unopy.scipy._hessian import _LagrangianHessianOperator
from unopy.scipy._objective import _Objective
from unopy.scipy._options import _set_options

if TYPE_CHECKING:
    from collections.abc import Iterable

AVAILABLE_METHODS = ["filterslp", "filtersqp", "funnelsqp", "ipopt"]
_FEASIBLE_SOLUTION_STATUSES = (
    unopy.SolutionStatus.FEASIBLE_KKT_POINT,
    unopy.SolutionStatus.FEASIBLE_FJ_POINT,
    unopy.SolutionStatus.FEASIBLE_SMALL_STEP,
)
"""The solution statuses of a feasible solution, required for a successful run."""


def _get_quasi_newton_hessian_model(
    hess: Callable[..., Any] | HessianUpdateStrategy | str | None,
) -> str | None:
    """Return the Uno Hessian model corresponding to a Hessian update strategy.

    Args:
        hess: The Hessian of the objective function, as in ``minimize``.

    Returns:
        The Uno Hessian model ("LBFGS" or "LSR1") if ``hess`` is a ``BFGS``
        or ``SR1`` instance, ``None`` if ``hess`` is a callable or ``None``.

    Raises:
        NotImplementedError: If ``hess`` is a finite-difference scheme
            or another Hessian update strategy.
        ValueError: If ``hess`` is not valid.
    """
    if hess is None or callable(hess):
        return None
    if isinstance(hess, SR1):
        return "LSR1"
    if isinstance(hess, BFGS):
        return "LBFGS"
    if isinstance(hess, HessianUpdateStrategy):
        raise NotImplementedError(
            f"The Hessian update strategy {hess!r} is not supported; "
            "use BFGS() or SR1()."
        )
    if _is_finite_difference_method(hess):
        raise NotImplementedError(
            f"The finite-difference Hessian {hess!r} is not supported; "
            "use a callable, BFGS() (Uno's L-BFGS) or SR1() (Uno's L-SR1)."
        )
    raise ValueError(
        "The Hessian 'hess' must be a callable, None, BFGS() or SR1(), "
        f"got {hess!r}."
    )


def minimize(
    fun: Callable[..., Any],
    x0: np.ndarray,
    args: tuple[Any, ...] = (),
    method: str = "filtersqp",
    jac: Callable[..., Any] | str | bool | None = None,
    hess: Callable[..., Any] | HessianUpdateStrategy | None = None,
    bounds: Bounds | Iterable[tuple[float | None, float | None]] | None = None,
    constraints: dict[str, Any]
    | NonlinearConstraint
    | LinearConstraint
    | Iterable[dict[str, Any] | NonlinearConstraint | LinearConstraint] = (),
    tol: float | None = None,
    callback: Callable[..., Any] | None = None,
    options: dict[str, Any] | None = None,
) -> OptimizeResult:
    """A scipy.optimize.minimize like interface for Uno.

    Minimization of scalar function of one or more variables under bounds and general constraints.

    Parameters
    ----------
    fun : callable
        The objective function to be minimized::

            fun(x, *args) -> float

        where ``x`` is a 1-D array with shape (n,) and ``args``
        is a tuple of the fixed parameters needed to completely
        specify the function.
        If ``jac=True``, ``fun`` returns the objective and its gradient::

            fun(x, *args) -> (float, array_like of shape (n,))

        and a single call to ``fun`` gives both values at a point.
        An array of size 1 (e.g. of shape (1,)) is accepted;
        a ``ValueError`` is raised for an array of another size.

        Suppose the callable has signature ``f0(x, *my_args, **my_kwargs)``, where
        ``my_args`` and ``my_kwargs`` are required positional and keyword arguments.
        Rather than passing ``f0`` as the callable, wrap it to accept
        only ``x``; e.g., pass ``fun=lambda x: f0(x, *my_args, **my_kwargs)`` as the
        callable, where ``my_args`` (tuple) and ``my_kwargs`` (dict) have been
        gathered before invoking this function.
    x0 : ndarray, shape (n,)
        Initial guess. Array of real elements of size (n,),
        where ``n`` is the number of independent variables.
    args : tuple, optional
        Extra arguments passed to the objective function and its
        derivatives ("fun", "jac" functions).
    method : str
        The UNO preset.  Should be one of

        - 'filtersqp'
        - 'filterslp'
        - "funnelsqp"
        - "ipopt"

    jac : {callable, '2-point', '3-point', 'cs', bool}, optional
        Method for computing the gradient vector.
        If it is a callable, it must be a function that returns the gradient
        vector::

            jac(x, *args) -> array_like, shape (n,)

        where ``x`` is an array with shape (n,) and ``args`` is a tuple with
        the fixed parameters.
        If ``jac`` is ``True``, ``fun`` returns the objective and its gradient;
        the last output of ``fun`` is stored with its point,
        so that one call gives both the objective and the gradient at a point.
        If ``jac`` is one of the finite-difference schemes '2-point' (forward),
        '3-point' (central) or 'cs' (complex step, ``fun`` must accept a complex
        ``x``), the gradient is approximated with the relative steps of SciPy;
        ``None`` (the default) and ``False`` mean '2-point'.
        The finite differences ignore the bounds:
        ``fun`` may be evaluated slightly outside them.
        '2-point' reuses the value of the objective computed by Uno at the point.

    hess : {callable, BFGS, SR1}, optional
        Hessian of the objective function. Signature::

            hess(x, *args) -> {ndarray (n, n), sparse matrix, LinearOperator}

        If it is a callable, the constraint Hessians are also read
        from the "NonlinearConstraint.hess" attributes.
        The constraint Hessian callable has the signature of SciPy::

            hess(x, v) -> {ndarray (n, n), sparse matrix, LinearOperator}

        returning the sum of the Hessians of the components of the constraint
        weighted by ``v``, the vector of the multipliers of this constraint.
        In that case, every nonlinear constraint must be a "NonlinearConstraint"
        with a callable "hess" ("LinearConstraint" has a zero Hessian);
        otherwise a "ValueError" is raised before the optimization.
        The Lagrangian Hessian is never formed as a dense matrix:
        Uno uses Hessian-vector products, computed term by term,
        so that a "scipy.sparse" matrix or a "LinearOperator" is never densified.
        The objective Hessian is computed once per point,
        the Hessian of a constraint once per point and multipliers,
        and "nhev" is the number of calls to "hess".

        If ``hess`` is a "BFGS" or "SR1" instance, Uno approximates the Lagrangian
        Hessian with its limited-memory quasi-Newton models, respectively
        L-BFGS (option "hessian_model" set to "LBFGS")
        and L-SR1 (option "hessian_model" set to "LSR1");
        the settings of the instance are not used
        and the option "quasi_newton_memory_size" sets the memory.
        The "hess" attributes of the constraints are then ignored.
        The method 'filterslp' uses no Hessian and ignores ``hess``.
        Without ``hess``, Uno uses its L-BFGS model.
        The finite-difference schemes ('2-point', '3-point', 'cs')
        are not supported for ``hess`` and raise a "NotImplementedError":
        Uno would call the Hessian-vector product at each iteration of its
        subproblem solver, i.e. a finite difference of the gradient
        of the Lagrangian per product, which is costly and inaccurate;
        use a quasi-Newton model instead.

    bounds : sequence or "Bounds", optional
        Bounds on variables. There are two ways to specify the bounds:

        1. Instance of "Bounds" class.
        2. Sequence of ``(min, max)`` pairs for each element in "x". None
           is used to specify no bound.

        The bounds are not modified.

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
        jac : {callable, '2-point', '3-point', 'cs'}, optional
            The Jacobian of "fun"; by default, '2-point'.
        args : sequence, optional
            Extra arguments to be passed to the function and Jacobian.

        Equality constraint means that the constraint function result is to
        be zero whereas inequality means that it is to be non-negative.

        The Jacobian of a "NonlinearConstraint" or a dictionary is either a callable
        or a finite-difference scheme ('2-point', '3-point' or 'cs'), as for ``jac``;
        a finite-difference Jacobian is dense
        and '2-point' reuses the value of the constraint computed by Uno at the point.
        A non-callable "NonlinearConstraint.hess" (e.g. the default "BFGS()")
        is ignored. The bounds "lb" and "ub" of a "NonlinearConstraint"
        are broadcast to the dimension of the constraint.
        When this dimension cannot be deduced from the bounds
        (scalar bounds or dictionary), the constraint function is evaluated at "x0";
        this value is reused for the first evaluation at "x0" by Uno.

        The Jacobian of a constraint and the matrix "A" of a "LinearConstraint"
        can be either array-like or a "scipy.sparse" matrix of shape (m, n),
        where "m" is the dimension of the constraint.
        A sparse Jacobian is kept sparse: only its stored entries are passed to Uno.
        The sparsity pattern is fixed for the whole optimization:
        it is the one of "A" or of the Jacobian at "x0"
        (all the entries for a dense Jacobian, the stored entries for a sparse one).
        To this end, the Jacobian of a nonlinear constraint given as a callable
        is evaluated at "x0" before the optimization; this value is reused for the first evaluation at "x0"
        by Uno. An entry of the pattern missing from a later evaluation is zero
        and a nonzero entry outside the pattern raises a "ValueError".

    tol : float, optional
        Tolerance for termination. When "tol" is specified, the selected
        minimization algorithm sets some relevant solver-specific tolerance(s)
        equal to "tol".

    callback : callable, optional
        A callable called by Uno at each accepted iterate, with the signature::

            callback(intermediate_result: OptimizeResult)

        where ``intermediate_result`` is a keyword parameter containing an
        "OptimizeResult" with the attributes

        - ``x``, the iterate,
        - ``fun``, the value of the objective at ``x``,
        - ``nit``, the number of the iterate, starting at 1,
        - ``constraint_multipliers``, the multipliers of the constraints
          concatenated in the order of ``constraints``,
        - ``lower_bound_multipliers`` and ``upper_bound_multipliers``,
          the multipliers of the bounds,
        - ``objective_multiplier``, the objective multiplier
          (0 in the feasibility restoration phase),
        - ``primal_feasibility``, ``stationarity`` and ``complementarity``,
          the residuals computed by Uno.

        The multipliers are the current estimates of Uno, with the sign convention
        of the result (see ``v`` below), possibly in the units of the functions
        scaled by Uno. ``fun`` is the value of the objective computed by Uno
        at ``x`` (stored with its point), so that no additional call to ``fun``
        is made in practice; if the last value was computed at another point,
        ``fun`` is evaluated at ``x``.
        Note that the name of the parameter must be ``intermediate_result``
        for the callback to be passed an "OptimizeResult";
        otherwise, the callback is called with a copy of ``x``::

            callback(xk)

        If the callback raises ``StopIteration``, Uno stops,
        the result is returned with ``success`` set to False
        and the message "`callback` raised `StopIteration`.";
        this exception is not re-raised.
        Any other exception raised by the callback stops Uno
        and is re-raised, as for the user functions.

    options : dict, optional
        A dictionary of solver options.

        The available options and their types are listed below.
        An int value is accepted for a float option.

        "primal_tolerance" : float
        "dual_tolerance" : float
        "loose_primal_tolerance" : float
        "loose_dual_tolerance" : float
        "loose_tolerance_iteration_threshold" : int
        "max_iterations" (or its aliases "maxiter" and "max_iter") : int
        "time_limit" : float
        "print_solution" : bool
        "unbounded_objective_threshold" : float
        "logger" : str
        "constraint_relaxation_strategy" : str
        "inequality_handling_method" : str
        "globalization_mechanism" : str
        "globalization_strategy" : str
        "hessian_model" : str
        "inertia_correction_strategy" : str
        "progress_norm" : str
        "residual_norm" : str
        "residual_scaling_threshold" : float
        "protect_actual_reduction_against_roundoff" : bool
        "protected_actual_reduction_macheps_coefficient" : float
        "print_subproblem" : bool
        "armijo_decrease_fraction" : float
        "armijo_tolerance" : float
        "switching_delta" : float
        "switching_infeasibility_exponent" : float
        "sufficient_infeasibility_decrease_ratio" : float
        "filter_type" : str
        "filter_beta" : float
        "filter_gamma" : float
        "filter_ubd" : float
        "filter_fact" : float
        "filter_capacity" : int
        "filter_sufficient_infeasibility_decrease_factor" : float
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
        "quasi_newton_memory_size" : int
        "quasi_newton_delta_lower_bound" : float
        "quasi_newton_delta_upper_bound" : float
        "LBFGS_max_skips_before_reset" : int
        "LSR1_pivot_max_magnitude" : float
        "regularization_failure_threshold" : float
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
        "barrier_function" : str
        "barrier_initial_parameter" : float
        "barrier_default_multiplier" : float
        "barrier_tau_min" : float
        "barrier_k_sigma" : float
        "barrier_k_mu" : float
        "barrier_theta_mu" : float
        "barrier_k_epsilon" : float
        "barrier_update_fraction" : float
        "barrier_regularization_exponent" : float
        "barrier_small_direction_factor" : float
        "barrier_push_variable_to_interior_k1" : float
        "barrier_push_variable_to_interior_k2" : float
        "barrier_damping_factor" : float
        "barrier_small_infeasibility_factor" : float
        "least_square_multiplier_max_norm" : float
        "BQPD_kmax_heuristic" : str
        "QP_solver" : str
        "LP_solver" : str
        "linear_solver" : str
        "option_file" : str
        "write_solution_to_file" : bool
        "SOC_max_iterations" : int
        "SOC_infeasibility_fraction" : float


    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        "OptimizeResult" for a description of other attributes.

        ``success`` is True if and only if Uno returns the optimization status
        "SUCCESS" with a feasible solution, i.e. the solution status
        "FEASIBLE_KKT_POINT", "FEASIBLE_FJ_POINT" or "FEASIBLE_SMALL_STEP";
        e.g. it is False at an infeasible stationary point
        or when the objective is unbounded.
        ``message`` is a string giving both the optimization status
        and the solution status, and ``status`` is the solution status.
        ``nfev`` is the number of calls to "fun",
        including the finite-difference ones.
        ``nhev`` is the number of calls to "hess" when it is a callable.

        The Lagrange multipliers follow the convention of SciPy's 'trust-constr':
        the Lagrangian is
        ``L(x, v, z) = f(x) + sum_i v_i^T c_i(x) + z^T x``
        and its gradient is zero at a KKT point,
        where ``c_i`` is the function of the i-th constraint
        (``A @ x`` for a "LinearConstraint"):

        - ``v`` is the list of the multipliers of the constraints,
          with one array per constraint in the order of ``constraints``,
          followed by the multipliers ``z`` of the bounds if ``bounds`` is given.
          A positive multiplier means that the upper bound is active
          and a negative one that the lower bound is active;
          a dictionary 'ineq' constraint ``g(x) >= 0`` has a nonpositive multiplier.
        - ``lower_bound_multipliers`` (nonpositive)
          and ``upper_bound_multipliers`` (nonnegative)
          are the multipliers of the lower and upper bounds of the variables,
          whose sum is ``z``.
        - ``constr`` is the list of the values of the constraint functions
          at the solution, with one array per constraint.

    Raises
    ------
    TypeError
        If the Jacobian of a constraint is neither a callable
        nor a finite-difference scheme.
    ValueError
        If "jac" is neither a callable, a boolean, None
        nor a finite-difference scheme,
        if "hess" is not valid, if the bounds are ill-shaped or inconsistent,
        if "hess" is a callable and a nonlinear constraint has no callable Hessian,
        or if a constraint Jacobian at "x0" has a wrong shape.
    NotImplementedError
        If "hess" is a finite-difference scheme
        or a Hessian update strategy other than "BFGS" and "SR1".
    BaseException
        The first exception raised by "fun", "jac", "hess", a constraint
        function or "callback" (except ``StopIteration``) during the optimization,
        including during the finite differences.
        Uno is stopped as soon as possible
        and the user functions are no longer called.

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

    # Step 1: Input validation
    x0 = np.atleast_1d(np.asarray(x0, dtype=float))
    n = len(x0)
    method_lower = method.lower()
    if method_lower not in AVAILABLE_METHODS:
        raise ValueError(
            f"Unknown method '{method}'. Must be one of {AVAILABLE_METHODS}"
        )
    if options is None:
        options = {}
    hessian_model = _get_quasi_newton_hessian_model(hess)

    # Step 2: Process bounds
    var_lb, var_ub = _get_variable_bounds(bounds, n)

    # Step 3: Build the objective
    objective = _Objective(fun, jac, args)

    # Step 4: Normalize and merge constraints
    normalized = _normalize_constraints(constraints, x0, hess)
    merged_constraints = _MergedConstraints(normalized)

    # Step 5: Build unopy.Model
    model = unopy.Model(
        unopy.PROBLEM_NONLINEAR, n, var_lb, var_ub, unopy.ZERO_BASED_INDEXING
    )
    guard = _EvaluationErrorGuard()
    model.set_objective(
        unopy.MINIMIZE,
        guard.wrap_objective(objective.evaluate),
        guard.wrap_inplace(objective.compute_gradient),
    )
    if merged_constraints.size > 0:
        model.set_constraints(
            merged_constraints.size,
            guard.wrap_inplace(merged_constraints.evaluate),
            merged_constraints.lower_bounds,
            merged_constraints.upper_bounds,
            merged_constraints.nnz,
            merged_constraints.row_indices,
            merged_constraints.col_indices,
            guard.wrap_inplace(merged_constraints.evaluate_jacobian),
        )
    model.set_initial_primal_iterate(list(x0))
    # The Lagrangian is f + v^T c, as in SciPy's trust-constr,
    # for both the Hessian operator and the multipliers of the result.
    model.set_lagrangian_sign_convention(unopy.MULTIPLIER_POSITIVE)

    # Step 5b: Lagrangian Hessian operator (if Hessian provided)
    hessian_operator = None
    if callable(hess):
        hessian_operator = _LagrangianHessianOperator(hess, args, normalized)
        model.set_lagrangian_hessian_operator(
            guard.wrap_inplace(hessian_operator, fill_with_zeros=False)
        )

    # Step 6: Configure and run solver
    uno_solver = unopy.UnoSolver()
    uno_solver.set_preset(method_lower)
    if hessian_model is not None and method_lower != "filterslp":
        # Set before the user options, which can override it.
        # The SLP method uses no Hessian.
        uno_solver.set_option("hessian_model", hessian_model)
    _set_options(uno_solver, tol, options)
    iteration_callback = None
    if callback is None:
        uno_solver.set_termination_callback(guard.terminate)
    else:
        iteration_callback = _IterationCallback(callback, n, objective, guard)
        uno_solver.set_notify_acceptable_iterate_callback(iteration_callback.notify)
        uno_solver.set_termination_callback(iteration_callback.terminate)

    result = uno_solver.optimize(model)
    guard.raise_exception()

    # Step 7: Return OptimizeResult
    # Uno may return the status SUCCESS at an infeasible stationary point.
    success = (
        result.optimization_status == unopy.OptimizationStatus.SUCCESS
        and result.solution_status in _FEASIBLE_SOLUTION_STATUSES
    )
    message = (
        f"Optimization status: {result.optimization_status.name}; "
        f"solution status: {result.solution_status.name}."
    )
    if iteration_callback is not None and iteration_callback.stopped:
        success = False
        message = STOP_ITERATION_MESSAGE
    if hessian_operator is None:
        nhev = result.number_hessian_evaluations
    else:
        nhev = hessian_operator.evaluations

    # Split the merged multipliers and constraint values per user constraint.
    constraint_multipliers = np.asarray(result.constraint_dual_solution, dtype=float)
    constraint_values = np.asarray(result.solution_constraints, dtype=float)
    lower_bound_multipliers = np.asarray(result.lower_bound_dual_solution, dtype=float)
    upper_bound_multipliers = np.asarray(result.upper_bound_dual_solution, dtype=float)
    offsets = np.cumsum([0] + [item.size for item in normalized])
    v = [
        constraint_multipliers[start:stop].copy()
        for start, stop in zip(offsets[:-1], offsets[1:])
    ]
    constr = [
        constraint_values[start:stop].copy()
        for start, stop in zip(offsets[:-1], offsets[1:])
    ]
    if bounds is not None:
        v.append(lower_bound_multipliers + upper_bound_multipliers)

    return OptimizeResult(
        x=result.primal_solution,
        success=success,
        status=result.solution_status,
        message=message,
        fun=result.solution_objective,
        nfev=objective.evaluations,
        njev=result.number_jacobian_evaluations,
        nhev=nhev,
        nit=result.number_iterations,
        maxcv=result.solution_primal_feasibility,
        v=v,
        constr=constr,
        lower_bound_multipliers=lower_bound_multipliers,
        upper_bound_multipliers=upper_bound_multipliers,
    )
