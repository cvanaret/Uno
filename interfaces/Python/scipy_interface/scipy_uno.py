# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
import numpy as np
import scipy.sparse as sp
import unopy
from numbers import Number
from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint, OptimizeResult
from typing import Iterable, Sequence, Any

Inf = float("inf")
AVAILABLE_METHODS = ["filterslp", "filtersqp", "funnelsqp", "ipopt"]
_MAX_ITERATIONS_ALIASES = ("maxiter", "max_iter")
"""The aliases of the Uno option "max_iterations"."""
_FEASIBLE_SOLUTION_STATUSES = (
    unopy.SolutionStatus.FEASIBLE_KKT_POINT,
    unopy.SolutionStatus.FEASIBLE_FJ_POINT,
    unopy.SolutionStatus.FEASIBLE_SMALL_STEP,
)
"""The solution statuses of a feasible solution, required for a successful run."""


def _process_bounds(bounds, n):
    """Return the lower and upper bounds of the variables as new float arrays.

    Args:
        bounds: Either a "Bounds" instance or a sequence of "(min, max)" pairs
            with "None" meaning no bound. It is not modified.
        n: The number of variables.

    Returns:
        The lower and upper bounds, each of shape (n,).

    Raises:
        ValueError: If the bounds are ill-shaped or a lower bound is greater
            than the corresponding upper bound.
    """
    if isinstance(bounds, Bounds):
        lb = np.asarray(bounds.lb, dtype=float)
        ub = np.asarray(bounds.ub, dtype=float)
        try:
            lb = np.array(np.broadcast_to(lb, (n,)), dtype=float)
            ub = np.array(np.broadcast_to(ub, (n,)), dtype=float)
        except ValueError:
            raise ValueError(
                f"The bounds must be broadcastable to the shape ({n},) of x0, "
                f"got lb of shape {lb.shape} and ub of shape {ub.shape}."
            ) from None
    else:
        pairs = list(bounds)
        if len(pairs) != n or any(np.ndim(pair) != 1 or len(pair) != 2 for pair in pairs):
            raise ValueError(
                f"The bounds must be a Bounds instance or a sequence of {n} "
                "(min, max) pairs."
            )
        lb = np.array([-Inf if low is None else low for low, _ in pairs], dtype=float)
        ub = np.array([Inf if high is None else high for _, high in pairs], dtype=float)

    if np.any(lb > ub):
        raise ValueError(
            "A lower bound is greater than the corresponding upper bound."
        )
    return lb, ub


def _reuse_first_value(function, x0, value):
    """Wrap a function to reuse its value at x0 for its first call at x0.

    Args:
        function: The function ``function(x) -> value``.
        x0: The point at which the function has already been evaluated.
        value: The value of the function at x0.

    Returns:
        The function returning ``value`` if its first call is at ``x0``.
    """
    cache = [(x0, value)]

    def wrapped_function(x):
        if cache:
            x_cached, value_cached = cache.pop()
            if np.array_equal(x, x_cached):
                return value_cached
        return function(x)

    return wrapped_function


def _set_option(uno_solver, name, value):
    """Set an option of Uno.

    The bindings choose the type of the option from the Python type of the value,
    so an integer value of a floating-point option is converted to a float.

    Args:
        uno_solver: The Uno solver.
        name: The name of the option.
        value: The value of the option.
    """
    try:
        uno_solver.set_option(name, value)
    except TypeError as error:
        if isinstance(value, (bool, np.bool_)) or not isinstance(
            value, (int, np.integer)
        ):
            raise
        try:
            uno_solver.set_option(name, float(value))
        except TypeError:
            raise error from None


class _JacobianBlock:
    """The Jacobian of a block of constraints with a sparsity pattern fixed for the run.

    The sparsity pattern is deduced from a first value of the Jacobian:
    all the entries of a dense Jacobian
    or the stored entries of a ``scipy.sparse`` Jacobian.
    The nonzeros are ordered row by row and, within a row, by increasing column.
    A sparse Jacobian is never converted to a dense array.
    """

    def __init__(self, value, x0, m, index, jacobian=None):
        """
        Args:
            value: The value of the Jacobian at ``x0``,
                either array-like or a ``scipy.sparse`` matrix.
            x0: The point at which ``value`` has been computed.
            m: The number of constraints of the block.
            index: The index of the constraint in the user's constraints,
                used in the error messages.
            jacobian: The function ``jacobian(x) -> value`` computing the Jacobian.
                If ``None``, the Jacobian is constant and equal to ``value``.
        """
        n = len(x0)
        self._index = index
        self._shape = (m, n)
        self._is_sparse = sp.issparse(value)
        if self._is_sparse:
            value = self._to_canonical_csr(value, copy=True)
            self.rows = np.repeat(np.arange(m), np.diff(value.indptr))
            self.cols = value.indices.astype(np.int64)
            # The linear indices of the nonzeros, sorted in increasing order.
            self._keys = self.rows * n + self.cols
        else:
            value = np.array(self._check_dense_size(value), dtype=float)
            self.rows = np.repeat(np.arange(m), n)
            self.cols = np.tile(np.arange(n), m)
        self.nnz = len(self.rows)
        """The number of nonzeros of the sparsity pattern."""

        if jacobian is None:
            self._constant_values = np.empty(self.nnz)
            self._write_values(value, self._constant_values)
            self._jacobian = None
        else:
            self._constant_values = None
            # The value at x0 is reused for the first evaluation by Uno.
            self._jacobian = _reuse_first_value(jacobian, x0, value)

    def _check_dense_size(self, value):
        """Return a dense Jacobian as an array after checking its size.

        Args:
            value: The dense Jacobian.

        Returns:
            The Jacobian as an array of floats.

        Raises:
            ValueError: If the Jacobian does not have the expected size.
        """
        value = np.asarray(value, dtype=float)
        if value.size != self._shape[0] * self._shape[1]:
            raise ValueError(
                f"The Jacobian of the constraint {self._index} must have the shape "
                f"{self._shape}, got an array of shape {value.shape}."
            )
        return value

    def _to_canonical_csr(self, value, copy=False):
        """Convert a sparse Jacobian to a CSR matrix in canonical format.

        Args:
            value: The sparse Jacobian.
            copy: Whether to copy the data.

        Returns:
            The CSR matrix with sorted indices and no duplicates.

        Raises:
            ValueError: If the Jacobian does not have the expected shape.
        """
        if value.shape != self._shape:
            raise ValueError(
                f"The Jacobian of the constraint {self._index} must have the shape "
                f"{self._shape}, got a sparse matrix of shape {value.shape}."
            )
        value = sp.csr_matrix(value, dtype=float, copy=copy)
        if not value.has_canonical_format:
            if not copy:
                value = value.copy()
            value.sum_duplicates()
        return value

    def _write_values(self, value, values):
        """Write the values of the Jacobian in the order of the sparsity pattern.

        Args:
            value: The value of the Jacobian.
            values: The output vector of size ``nnz``.

        Raises:
            ValueError: If the Jacobian has a nonzero entry
                outside the sparsity pattern.
        """
        if not self._is_sparse:
            if sp.issparse(value):
                value = self._to_canonical_csr(value).toarray()
            values[:] = self._check_dense_size(value).ravel()
            return

        n = self._shape[1]
        if sp.issparse(value):
            value = self._to_canonical_csr(value)
            keys = (
                np.repeat(np.arange(self._shape[0]), np.diff(value.indptr)) * n
                + value.indices
            )
            data = value.data
        else:
            value = self._check_dense_size(value).ravel()
            keys = np.flatnonzero(value)
            data = value[keys]

        positions = np.searchsorted(self._keys, keys)
        in_pattern = positions < self.nnz
        in_pattern[in_pattern] = self._keys[positions[in_pattern]] == keys[in_pattern]
        outside = ~in_pattern & (data != 0.0)
        if np.any(outside):
            row, col = divmod(int(keys[np.argmax(outside)]), n)
            raise ValueError(
                f"The Jacobian of the constraint {self._index} has a nonzero entry "
                f"at ({row}, {col}) outside its sparsity pattern, "
                "which is the one of its Jacobian at x0."
            )
        values[:] = 0.0
        values[positions[in_pattern]] = data[in_pattern]

    def write(self, x, values):
        """Evaluate the Jacobian and write its values in the order of the pattern.

        Args:
            x: The point at which the Jacobian is evaluated.
            values: The output vector of size ``nnz``.
        """
        if self._jacobian is None:
            values[:] = self._constant_values
        else:
            self._write_values(self._jacobian(x), values)


class _EvaluationErrorGuard:
    """Catch the exceptions raised by the user functions during an optimization.

    Uno converts any exception raised in a Python callback into an evaluation error
    and keeps iterating (e.g. by backtracking), so the original exception would be lost.
    Instead, the first exception raised by a user function is stored,
    the callback returns zeros
    (except the Hessian operator, which raises the exception to Uno),
    and the termination callback asks Uno to stop.
    Once an exception has been stored, the user functions are no longer called
    and the callbacks keep returning zeros (or raising for the Hessian operator).
    Zeros are used rather than NaN values
    because NaN values can make the LP/QP solvers loop forever
    before the termination callback is called;
    the values returned after a failure are discarded anyway
    since the exception is re-raised once Uno has returned.
    """

    def __init__(self):
        self.exception = None
        """The first exception raised by a user function, if any."""

    def wrap_objective(self, callback):
        """Wrap a callback returning the objective value.

        Args:
            callback: The callback ``callback(x) -> float``.

        Returns:
            The callback returning 0.0 instead of raising an exception.
        """

        def guarded_callback(x):
            if self.exception is None:
                try:
                    return callback(x)
                except BaseException as exception:
                    self.exception = exception
            return 0.0

        return guarded_callback

    def wrap_inplace(self, callback, fill_with_zeros=True):
        """Wrap a callback writing its result in its last argument.

        Args:
            callback: The callback ``callback(*inputs, output) -> None``.
            fill_with_zeros: Whether to fill the output with zeros after a failure.
                Otherwise, the exception is raised to Uno,
                which stops with an evaluation error.
                This is used for the Hessian operator.

        Returns:
            The callback filling the output with zeros instead of raising an exception.
        """

        def guarded_callback(*arguments):
            if self.exception is None:
                try:
                    callback(*arguments)
                    return
                except BaseException as exception:
                    self.exception = exception
            if not fill_with_zeros:
                raise self.exception
            arguments[-1][:] = 0.0

        return guarded_callback

    def terminate(self, *_):
        """Return whether Uno must stop because a user function raised an exception.

        Args:
            *_: The arguments passed by Uno (iterate, multipliers and residuals).

        Returns:
            Whether a user function raised an exception.
        """
        return self.exception is not None

    def raise_exception(self):
        """Re-raise the exception raised by a user function, if any.

        Raises:
            BaseException: The exception raised by a user function.
        """
        if self.exception is not None:
            raise self.exception


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
        An array of size 1 (e.g. of shape (1,)) is accepted;
        a ""ValueError"" is raised for an array of another size.

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

    jac : {callable}
        Method for computing the gradient vector.
        It must be a function that returns the gradient
        vector::

            jac(x, *args) -> array_like, shape (n,)

        where ""x"" is an array with shape (n,) and ""args"" is a tuple with
        the fixed parameters.
        "jac=True" (i.e. "fun" returning the objective and its gradient)
        and finite-difference schemes are not supported.

    hess : callable, optional
        Hessian of the objective function. Signature: hess(x, *args) -> ndarray (n, n).
        If provided, constraint Hessians are also read from NonlinearConstraint.hess
        attributes. The constraint Hessian callable has signature:
        hess(x, v) -> ndarray (n, n), where v is the vector of constraint multipliers.
        In that case, every nonlinear constraint must be a "NonlinearConstraint"
        with a callable "hess" ("LinearConstraint" has a zero Hessian);
        otherwise a "ValueError" is raised before the optimization.
        The Lagrangian Hessian is computed once per point
        and "nhev" is the number of calls to "hess".

    bounds : sequence or "Bounds", optional
        Bounds on variables. There are two ways to specify the bounds:

        1. Instance of "Bounds" class.
        2. Sequence of ""(min, max)"" pairs for each element in "x". None
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
        jac : callable
            The Jacobian of "fun" (required).
        args : sequence, optional
            Extra arguments to be passed to the function and Jacobian.

        Equality constraint means that the constraint function result is to
        be zero whereas inequality means that it is to be non-negative.

        The Jacobian of a "NonlinearConstraint" must be a callable:
        finite-difference schemes such as '2-point' are not supported.
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
        To this end, the Jacobian of a nonlinear constraint is evaluated at "x0"
        before the optimization; this value is reused for the first evaluation at "x0"
        by Uno. An entry of the pattern missing from a later evaluation is zero
        and a nonzero entry outside the pattern raises a "ValueError".

    tol : float, optional
        Tolerance for termination. When "tol" is specified, the selected
        minimization algorithm sets some relevant solver-specific tolerance(s)
        equal to "tol".

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
        The optimization result represented as a ""OptimizeResult"" object.
        Important attributes are: ""x"" the solution array, ""success"" a
        Boolean flag indicating if the optimizer exited successfully and
        ""message"" which describes the cause of the termination. See
        "OptimizeResult" for a description of other attributes.

        ""success"" is True if and only if Uno returns the optimization status
        "SUCCESS" with a feasible solution, i.e. the solution status
        "FEASIBLE_KKT_POINT", "FEASIBLE_FJ_POINT" or "FEASIBLE_SMALL_STEP";
        e.g. it is False at an infeasible stationary point
        or when the objective is unbounded.
        ""message"" is a string giving both the optimization status
        and the solution status, and ""status"" is the solution status.
        ""nhev"" is the number of calls to "hess" when it is provided.

    Raises
    ------
    TypeError
        If the Jacobian of a constraint is not a callable.
    ValueError
        If "jac" is not a callable, if the bounds are ill-shaped or inconsistent,
        if "hess" is provided and a nonlinear constraint has no callable Hessian,
        or if a constraint Jacobian at "x0" has a wrong shape.
    BaseException
        The first exception raised by "fun", "jac", "hess" or a constraint
        function during the optimization. Uno is stopped as soon as possible
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

    # Step 2: Process bounds
    if bounds is not None:
        var_lb, var_ub = _process_bounds(bounds, n)
    else:
        var_lb = [-Inf] * n
        var_ub = [Inf] * n

    # Step 3: Build objective callbacks
    def objective_callback(x):
        value = np.asarray(fun(x, *args), dtype=float)
        if value.size != 1:
            raise ValueError(
                "The objective function must return a scalar, "
                f"got an array of shape {value.shape}."
            )
        return value.item()

    if jac is None:
        raise ValueError("Objective gradient is not defined.")
    if not callable(jac):
        raise ValueError(
            f"The objective gradient 'jac' must be a callable, got {jac!r}; "
            "'jac=True' (fun returning the objective and its gradient) "
            "and finite differences are not supported."
        )

    def gradient_callback(x, gradient):
        gradient[:] = np.asarray(jac(x, *args), dtype=float).ravel()

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
            if not callable(c_jac):
                raise TypeError(
                    f"The Jacobian of the constraint {i} must be a callable, "
                    f"got {c_jac!r}; finite differences are not supported."
                )
            c_hess = getattr(con, "hess", None)
            if not callable(c_hess):
                # E.g. the default BFGS() strategy or a finite-difference scheme.
                c_hess = None
            # The size of the constraint is given by its bounds if they are not scalars;
            # otherwise, the constraint is evaluated at x0
            # and this value is reused for the first evaluation by Uno.
            bounds_size = np.broadcast(
                np.asarray(con.lb, dtype=float), np.asarray(con.ub, dtype=float)
            ).size
            if bounds_size > 1:
                m_i = bounds_size
            else:
                value = c_fun(x0.copy())
                m_i = len(np.atleast_1d(value))
                c_fun = _reuse_first_value(c_fun, x0, value)
            c_lb = np.broadcast_to(
                np.atleast_1d(np.asarray(con.lb, dtype=float)), (m_i,)
            ).copy()
            c_ub = np.broadcast_to(
                np.atleast_1d(np.asarray(con.ub, dtype=float)), (m_i,)
            ).copy()
            # The Jacobian at x0 gives the sparsity pattern;
            # it is reused for the first evaluation by Uno.
            c_jac = _JacobianBlock(c_jac(x0.copy()), x0, m_i, i, jacobian=c_jac)
            normalized.append((c_fun, c_jac, c_hess, c_lb, c_ub, m_i))
        elif isinstance(con, LinearConstraint):
            if sp.issparse(con.A):
                A = sp.csr_matrix(con.A, dtype=float)
            else:
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
            # A constant Jacobian whose sparsity pattern is the one of A.
            c_jac = _JacobianBlock(A, x0, m_i, i)
            c_hess = None  # Linear constraints have zero Hessian
            normalized.append((c_fun, c_jac, c_hess, c_lb, c_ub, m_i))
        elif isinstance(con, dict):
            c_fun_raw = con["fun"]
            c_jac_raw = con.get("jac")
            if not callable(c_jac_raw):
                raise TypeError(
                    f"The 'jac' of the constraint {i} must be a callable, "
                    f"got {c_jac_raw!r}."
                )
            c_args = con.get("args", ())
            c_type = con["type"]
            # Evaluate to get dimension; this value is reused for the first evaluation by Uno.
            val = np.atleast_1d(np.asarray(c_fun_raw(x0.copy(), *c_args), dtype=float))
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
            c_fun = _reuse_first_value(
                lambda x, _f=c_fun_raw, _a=c_args: np.atleast_1d(
                    np.asarray(_f(x, *_a), dtype=float)
                ),
                x0,
                val,
            )
            c_jac = lambda x, _j=c_jac_raw, _a=c_args: _j(x, *_a)
            # The Jacobian at x0 gives the sparsity pattern;
            # it is reused for the first evaluation by Uno.
            c_jac = _JacobianBlock(c_jac(x0.copy()), x0, m_i, i, jacobian=c_jac)
            c_hess = None  # Dict constraints don't support Hessians
            normalized.append((c_fun, c_jac, c_hess, c_lb, c_ub, m_i))
        else:
            raise TypeError(f"Unsupported constraint type: {type(con)}")

        if hess is not None and c_hess is None and not isinstance(
            con, LinearConstraint
        ):
            raise ValueError(
                f"The constraint {i} has no callable Hessian, "
                "which is required when 'hess' is provided; "
                "use a NonlinearConstraint with a callable 'hess'."
            )

    # 4c: Merge
    total_constraints = sum(item[5] for item in normalized)
    if total_constraints > 0:
        all_lb = np.concatenate([item[3] for item in normalized])
        all_ub = np.concatenate([item[4] for item in normalized])
    else:
        all_lb = np.array([])
        all_ub = np.array([])

    # 4d: Build merged callbacks
    def constraint_callback(x, constraint_values):
        offset = 0
        for c_fun_i, _, _, _, _, m_i in normalized:
            constraint_values[offset:offset + m_i] = np.asarray(c_fun_i(x), dtype=float).ravel()
            offset += m_i

    def constraint_jacobian_callback(x, jacobian_values):
        # The values of each block are written in the order of the sparsity pattern
        # below, without building the dense merged Jacobian.
        start = 0
        for _, block, _, _, _, _ in normalized:
            block.write(x, jacobian_values[start:start + block.nnz])
            start += block.nnz

    # 4e: Sparsity pattern (COO), concatenation of the patterns of the blocks
    if total_constraints > 0:
        nnz_jac = sum(item[1].nnz for item in normalized)
        row_offsets = np.cumsum([0] + [item[5] for item in normalized[:-1]])
        row_indices = np.concatenate(
            [item[1].rows + offset for item, offset in zip(normalized, row_offsets)]
        ).tolist()
        col_indices = np.concatenate([item[1].cols for item in normalized]).tolist()
    else:
        nnz_jac = 0
        row_indices = []
        col_indices = []

    # Step 5: Build unopy.Model
    model = unopy.Model(
        unopy.PROBLEM_NONLINEAR, n, var_lb, var_ub, unopy.ZERO_BASED_INDEXING
    )
    guard = _EvaluationErrorGuard()
    model.set_objective(
        unopy.MINIMIZE,
        guard.wrap_objective(objective_callback),
        guard.wrap_inplace(gradient_callback),
    )
    if total_constraints > 0:
        model.set_constraints(
            total_constraints,
            guard.wrap_inplace(constraint_callback),
            all_lb,
            all_ub,
            nnz_jac,
            row_indices,
            col_indices,
            guard.wrap_inplace(constraint_jacobian_callback),
        )
    model.set_initial_primal_iterate(list(x0))

    # Step 5b: Lagrangian Hessian operator (if Hessian provided)
    if hess is not None:
        # The last Lagrangian Hessian and the point (x, obj_mult, multipliers)
        # at which it was evaluated, reused by the successive Hessian-vector
        # products computed by the subproblem solver at the same point.
        hessian_cache = {}
        # The number of evaluations of the user Hessian:
        # Uno does not count the evaluations of a Hessian operator.
        hessian_evaluations = [0]

        def compute_lagrangian_hessian(x, obj_mult, multipliers):
            # Objective Hessian contribution
            hessian_evaluations[0] += 1
            H =float(obj_mult) * np.atleast_2d(np.asarray(hess(x, *args), dtype=float))

            # Constraint Hessian contributions
            offset = 0
            for _, _, c_hess_i, _, _, m_i in normalized:
                if c_hess_i is not None:
                    mult_slice = np.asarray(multipliers[offset:offset + m_i], dtype=float)
                    H = H + np.atleast_2d(np.asarray(c_hess_i(x, mult_slice), dtype=float))
                offset += m_i
            return H

        def hessian_operator_callback(x, evaluate_at_x, obj_mult,
                                      multipliers, vector, result):
            # Uno passes views on its own buffers: copy them to build the cache key.
            point = (
                np.array(x, dtype=float),
                float(obj_mult),
                np.array(multipliers, dtype=float),
            )
            cached_point = hessian_cache.get("point")
            if cached_point is None or (
                evaluate_at_x
                and not (
                    np.array_equal(point[0], cached_point[0])
                    and point[1] == cached_point[1]
                    and np.array_equal(point[2], cached_point[2])
                )
            ):
                # Store the point only once the Hessian has been computed,
                # so that a failed evaluation is not cached.
                hessian_cache["H"] = compute_lagrangian_hessian(*point)
                hessian_cache["point"] = point

            # Compute H @ vector
            result[:] = hessian_cache["H"] @ vector

        model.set_lagrangian_hessian_operator(
            guard.wrap_inplace(hessian_operator_callback, fill_with_zeros=False)
        )
        model.set_lagrangian_sign_convention(unopy.MULTIPLIER_POSITIVE)


    # Step 6: Configure and run solver
    uno_solver = unopy.UnoSolver()
    uno_solver.set_preset(method_lower)
    if tol is not None:
        uno_solver.set_option("primal_tolerance", float(tol))
        uno_solver.set_option("dual_tolerance", float(tol))
    for key, value in options.items():
        if key in _MAX_ITERATIONS_ALIASES:
            key = "max_iterations"
        _set_option(uno_solver, key, value)

    uno_solver.set_termination_callback(guard.terminate)

    result = uno_solver.optimize(model)
    guard.raise_exception()

    # Step 7: Return OptimizeResult
    # Uno may return the status SUCCESS at an infeasible stationary point.
    success = (
        result.optimization_status == unopy.OptimizationStatus.SUCCESS
        and result.solution_status in _FEASIBLE_SOLUTION_STATUSES
    )
    if hess is None:
        nhev = result.number_hessian_evaluations
    else:
        nhev = hessian_evaluations[0]
    return OptimizeResult(
        x=result.primal_solution,
        success=success,
        status=result.solution_status,
        message=(
            f"Optimization status: {result.optimization_status.name}; "
            f"solution status: {result.solution_status.name}."
        ),
        fun=result.solution_objective,
        nfev=result.number_objective_evaluations,
        njev=result.number_jacobian_evaluations,
        nhev=nhev,
        nit=result.number_iterations,
        maxcv=result.solution_primal_feasibility,
    )
