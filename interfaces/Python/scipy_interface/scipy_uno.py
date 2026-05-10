# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
import numpy as np
import unopy
from numbers import Number
from scipy.optimize import LinearConstraint, NonlinearConstraint, OptimizeResult
from scipy.optimize._minimize import standardize_bounds, _validate_bounds, Bounds
from typing import Iterable, Sequence, Any

Inf = float("inf")
AVAILABLE_METHODS = ["filterslp", "filtersqp", "funnelsqp", "ipopt"]


def minimize(
    fun: callable,
    x0: np.ndarray,
    args: tuple = (),
    method: str = "filtersqp",
    jac: callable = None,
    hess: callable = None,
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

    hess : callable, optional
        Hessian of the objective function. Signature: hess(x, *args) -> ndarray (n, n).
        If provided, constraint Hessians are also read from NonlinearConstraint.hess
        attributes. The constraint Hessian callable has signature:
        hess(x, v) -> ndarray (n, n), where v is the vector of constraint multipliers.

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

    # Step 2: Process bounds
    if bounds is not None:
        bounds = standardize_bounds(bounds, x0, "new")
        _validate_bounds(bounds, x0, method)
        var_lb = list(bounds.lb)
        var_ub = list(bounds.ub)
    else:
        var_lb = [-Inf] * n
        var_ub = [Inf] * n

    # Helper to convert unopy array to numpy
    def _to_numpy(x, size):
        return np.array([float(x[i]) for i in range(size)])

    # Step 3: Build objective callbacks
    def objective_callback(number_variables, x, objective_value, user_data):
        x_np = _to_numpy(x, number_variables)
        objective_value[0] = float(fun(x_np, *args))
        return 0

    if jac is not None:
        def gradient_callback(number_variables, x, gradient, user_data):
            x_np = _to_numpy(x, number_variables)
            g = np.asarray(jac(x_np, *args), dtype=float).ravel()
            for i in range(number_variables):
                gradient[i] = float(g[i])
            return 0
    else:
        raise ValueError("Objective gradient is not defined.")

    # Step 4: Normalize and merge constraints
    # 4a: Normalize input to a list
    if isinstance(constraints, dict) or isinstance(
        constraints, (NonlinearConstraint, LinearConstraint)
    ):
        constraints = [constraints]
    else:
        constraints = list(constraints)

    # 4b: Convert each constraint to (c_fun, c_jac, c_hess, c_lb, c_ub, m_i)
    normalized = []
    for i, con in enumerate(constraints):
        if isinstance(con, NonlinearConstraint):
            c_fun = con.fun
            c_jac = con.jac
            c_hess = getattr(con, "hess", None)
            c_lb = np.atleast_1d(np.asarray(con.lb, dtype=float))
            c_ub = np.atleast_1d(np.asarray(con.ub, dtype=float))
            m_i = len(np.atleast_1d(c_fun(x0)))
            normalized.append((c_fun, c_jac, c_hess, c_lb, c_ub, m_i))
        elif isinstance(con, LinearConstraint):
            A = np.atleast_2d(np.asarray(con.A, dtype=float))
            m_i = A.shape[0]
            c_lb = np.broadcast_to(
                np.atleast_1d(np.asarray(con.lb, dtype=float)), (m_i,)
            ).copy()
            c_ub = np.broadcast_to(
                np.atleast_1d(np.asarray(con.ub, dtype=float)), (m_i,)
            ).copy()
            # Use default-arg binding to capture A
            c_fun = lambda x, _A=A: _A @ x
            c_jac = lambda x, _A=A: _A
            c_hess = None  # Linear constraints have zero Hessian
            normalized.append((c_fun, c_jac, c_hess, c_lb, c_ub, m_i))
        elif isinstance(con, dict):
            c_fun_raw = con["fun"]
            c_jac_raw = con.get("jac")
            c_args = con.get("args", ())
            c_type = con["type"]
            # Evaluate to get dimension
            val = np.atleast_1d(np.asarray(c_fun_raw(x0, *c_args), dtype=float))
            m_i = len(val)
            if c_type == "eq":
                c_lb = np.zeros(m_i)
                c_ub = np.zeros(m_i)
            elif c_type == "ineq":
                c_lb = np.zeros(m_i)
                c_ub = np.full(m_i, Inf)
            else:
                raise ValueError(
                    f"Unknown constraint type '{c_type}'. Must be 'eq' or 'ineq'"
                )
            # Use default-arg binding
            c_fun = lambda x, _f=c_fun_raw, _a=c_args: np.atleast_1d(
                np.asarray(_f(x, *_a), dtype=float)
            )
            if c_jac_raw is not None:
                c_jac = lambda x, _j=c_jac_raw, _a=c_args: np.atleast_2d(
                    np.asarray(_j(x, *_a), dtype=float)
                )
            else:
                c_jac = None
            c_hess = None  # Dict constraints don't support Hessians
            normalized.append((c_fun, c_jac, c_hess, c_lb, c_ub, m_i))
        else:
            raise TypeError(f"Unsupported constraint type: {type(con)}")

    # 4c: Merge
    total_constraints = sum(item[5] for item in normalized)
    if total_constraints > 0:
        all_lb = np.concatenate([item[3] for item in normalized])
        all_ub = np.concatenate([item[4] for item in normalized])
    else:
        all_lb = np.array([])
        all_ub = np.array([])

    # 4d: Build merged callbacks
    def constraint_callback(
        number_variables, number_constraints, x, constraint_values, user_data
    ):
        x_np = _to_numpy(x, number_variables)
        offset = 0
        for c_fun_i, _, _, _, _, m_i in normalized:
            val = np.atleast_1d(np.asarray(c_fun_i(x_np), dtype=float))
            for j in range(m_i):
                constraint_values[offset + j] = float(val[j])
            offset += m_i
        return 0

    def constraint_jacobian_callback(
        number_variables, number_jacobian_nonzeros, x, jacobian_values, user_data
    ):
        x_np = _to_numpy(x, number_variables)
        offset = 0
        for c_fun_i, c_jac_i, _, _, _, m_i in normalized:
            if c_jac_i is not None:
                jac_block = np.atleast_2d(np.asarray(c_jac_i(x_np), dtype=float))
            else:
                raise ValueError("Constraint Jacobian is not provided.")
            # Column-major (Fortran) ordering to match unopy convention
            flat = jac_block.ravel(order="F")
            for j in range(len(flat)):
                jacobian_values[offset + j] = float(flat[j])
            offset += len(flat)
        return 0

    # 4e: Dense sparsity pattern (column-major ordering)
    if total_constraints > 0:
        nnz_jac = total_constraints * n
        # Column-major: iterate columns first, then rows within each column
        col_indices = np.repeat(np.arange(n), total_constraints).tolist()
        row_indices = np.tile(np.arange(total_constraints), n).tolist()
    else:
        nnz_jac = 0
        row_indices = []
        col_indices = []

    # Step 5: Build unopy.Model
    model = unopy.Model(
        unopy.PROBLEM_NONLINEAR, n, var_lb, var_ub, unopy.ZERO_BASED_INDEXING
    )
    model.set_objective(unopy.MINIMIZE, objective_callback, gradient_callback)
    if total_constraints > 0:
        model.set_constraints(
            total_constraints,
            constraint_callback,
            list(all_lb),
            list(all_ub),
            nnz_jac,
            row_indices,
            col_indices,
            constraint_jacobian_callback,
        )
    model.set_initial_primal_iterate(list(x0))

    # Step 5b: Lagrangian Hessian operator (if Hessian provided)
    if hess is not None:
        def hessian_operator_callback(n_vars, n_cons, x, evaluate_at_x, obj_mult,
                                      multipliers, vector, result, user_data):
            x_np = _to_numpy(x, n_vars)
            v_np = _to_numpy(vector, n_vars)

            # Objective Hessian contribution
            H = float(obj_mult) * np.atleast_2d(np.asarray(hess(x_np, *args), dtype=float))

            # Constraint Hessian contributions
            offset = 0
            for _, _, c_hess_i, _, _, m_i in normalized:
                if c_hess_i is not None:
                    mult_slice = np.array([float(multipliers[offset + j]) for j in range(m_i)])
                    H += np.atleast_2d(np.asarray(c_hess_i(x_np, mult_slice), dtype=float))
                offset += m_i

            # Compute H @ vector and write element-by-element
            Hv = H @ v_np
            for i in range(n_vars):
                result[i] = float(Hv[i])
            return 0

        model.set_lagrangian_hessian_operator(
            hessian_operator_callback, unopy.MULTIPLIER_POSITIVE
        )

    # Step 6: Configure and run solver
    uno_solver = unopy.UnoSolver()
    uno_solver.set_preset(method_lower)
    if tol is not None:
        uno_solver.set_option("primal_tolerance", tol)
        uno_solver.set_option("dual_tolerance", tol)
    for key, value in options.items():
        if key == "maxiter":
            key = "max_iterations"
        uno_solver.set_option(key, value)

    result = uno_solver.optimize(model)

    # Step 7: Return OptimizeResult
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
