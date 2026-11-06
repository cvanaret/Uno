# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""The normalization and the merging of the constraints."""

from __future__ import annotations

from typing import Any
from typing import Callable
from typing import NamedTuple
from typing import TYPE_CHECKING

import numpy as np
import scipy.sparse as sp
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint

from unopy.scipy._finite_differences import FINITE_DIFFERENCE_METHODS
from unopy.scipy._finite_differences import _approximate_jacobian
from unopy.scipy._finite_differences import _is_finite_difference_method

if TYPE_CHECKING:
    from collections.abc import Iterable

Inf = float("inf")


def _reuse_first_value(
    function: Callable[[np.ndarray], Any], x0: np.ndarray, value: Any
) -> Callable[[np.ndarray], Any]:
    """Wrap a function to reuse its value at x0 for its first call at x0.

    Args:
        function: The function ``function(x) -> value``.
        x0: The point at which the function has already been evaluated.
        value: The value of the function at x0.

    Returns:
        The function returning ``value`` if its first call is at ``x0``.
    """
    cache = [(x0, value)]

    def wrapped_function(x: np.ndarray) -> Any:
        if cache:
            x_cached, value_cached = cache.pop()
            if np.array_equal(x, x_cached):
                return value_cached
        return function(x)

    return wrapped_function


class _LastValue:
    """A function storing its last value with its point."""

    def __init__(self, function: Callable[[np.ndarray], Any]) -> None:
        """
        Args:
            function: The function ``function(x) -> array_like``.
        """
        self._function = function
        self._x = None
        self._value = None

    def __call__(self, x: np.ndarray) -> Any:
        """Evaluate the function and store its value.

        Args:
            x: The point.

        Returns:
            The value of the function.
        """
        value = self._function(x)
        # Uno passes a view on its own buffer: copy it.
        self._x = np.array(x, dtype=float)
        self._value = np.array(value, dtype=float)
        return value

    def value_at(self, x: np.ndarray) -> np.ndarray:
        """Return the value of the function, reusing the last one if possible.

        Args:
            x: The point.

        Returns:
            The value of the function.
        """
        if self._x is None or not np.array_equal(x, self._x):
            self(x)
        return self._value


def _create_finite_difference_jacobian(
    function: Callable[[np.ndarray], Any],
    last_value: _LastValue,
    method: str,
) -> Callable[[np.ndarray], np.ndarray]:
    """Create a function approximating a Jacobian with finite differences.

    Args:
        function: The function ``function(x) -> array_like``
            used for the finite differences;
            with "cs", it must accept a complex ``x``.
        last_value: The function storing the last value computed for Uno,
            reused by "2-point".
        method: The finite-difference scheme, either "2-point", "3-point" or "cs".

    Returns:
        The function ``jacobian(x) -> ndarray (m, n)``.
    """

    def jacobian(x: np.ndarray) -> np.ndarray:
        value = last_value.value_at(x) if method == "2-point" else None
        return _approximate_jacobian(function, x, method, value)

    return jacobian


class _JacobianBlock:
    """The Jacobian of a block of constraints with a sparsity pattern fixed for the run.

    The sparsity pattern is deduced from a first value of the Jacobian:
    all the entries of a dense Jacobian
    or the stored entries of a ``scipy.sparse`` Jacobian.
    The nonzeros are ordered row by row and, within a row, by increasing column.
    A sparse Jacobian is never converted to a dense array.
    Without a first value, the sparsity pattern is dense.
    """

    def __init__(
        self,
        value: Any | None,
        x0: np.ndarray,
        m: int,
        index: int,
        jacobian: Callable[[np.ndarray], Any] | None = None,
    ) -> None:
        """
        Args:
            value: The value of the Jacobian at ``x0``,
                either array-like or a ``scipy.sparse`` matrix.
                If ``None``, the sparsity pattern is dense
                and ``jacobian`` must be given.
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
        if value is None:
            self._is_sparse = False
            self.rows = np.repeat(np.arange(m), n)
            self.cols = np.tile(np.arange(n), m)
            self.nnz = m * n
            self._constant_values = None
            self._jacobian = jacobian
            return

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

    def _check_dense_size(self, value: Any) -> np.ndarray:
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

    def _to_canonical_csr(self, value: Any, copy: bool = False) -> sp.csr_matrix:
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

    def _write_values(self, value: Any, values: np.ndarray) -> None:
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

    def write(self, x: np.ndarray, values: np.ndarray) -> None:
        """Evaluate the Jacobian and write its values in the order of the pattern.

        Args:
            x: The point at which the Jacobian is evaluated.
            values: The output vector of size ``nnz``.
        """
        if self._jacobian is None:
            values[:] = self._constant_values
        else:
            self._write_values(self._jacobian(x), values)


class _Constraint(NamedTuple):
    """A block of constraints normalized for Uno."""

    fun: Callable[[np.ndarray], Any]
    """The function ``fun(x) -> value`` computing the constraint."""

    jacobian: _JacobianBlock
    """The Jacobian of the constraint."""

    hessian: Callable[[np.ndarray, np.ndarray], Any] | None
    """The function ``hessian(x, v)``, if any,
    returning ``sum_i v_i * Hessian_i(x)``
    where ``v`` is the vector of constraint multipliers,
    as a dense array, a ``scipy.sparse`` matrix or a ``LinearOperator``."""

    lower_bound: np.ndarray
    """The lower bound of the constraint."""

    upper_bound: np.ndarray
    """The upper bound of the constraint."""

    size: int
    """The dimension of the constraint."""


def _normalize_constraints(
    constraints: dict[str, Any]
    | NonlinearConstraint
    | LinearConstraint
    | Iterable[dict[str, Any] | NonlinearConstraint | LinearConstraint],
    x0: np.ndarray,
    hess: Callable[..., Any] | None,
) -> list[_Constraint]:
    """Normalize the constraints.

    Args:
        constraints: The constraints, as in ``minimize``.
        x0: The initial guess.
        hess: The Hessian of the objective function, if any;
            the Hessians of the nonlinear constraints are required
            only if it is a callable.

    Returns:
        The normalized constraints.

    Raises:
        TypeError: If the Jacobian of a constraint is neither a callable
            nor a finite-difference scheme, or a constraint has an unsupported type.
        ValueError: If a dictionary has an unknown constraint type,
            if ``hess`` is a callable and a nonlinear constraint has no callable Hessian,
            or if a constraint Jacobian at ``x0`` has a wrong shape.
    """
    # Normalize input to a list
    if isinstance(constraints, dict) or isinstance(
        constraints, (NonlinearConstraint, LinearConstraint)
    ):
        constraints = [constraints]
    else:
        constraints = list(constraints)

    # Convert each constraint to (c_fun, c_jac, c_hess, c_lb, c_ub, m_i)
    normalized = []
    for i, con in enumerate(constraints):
        if isinstance(con, NonlinearConstraint):
            c_fun = con.fun
            c_jac = con.jac
            if not (callable(c_jac) or _is_finite_difference_method(c_jac)):
                raise TypeError(
                    f"The Jacobian of the constraint {i} must be a callable "
                    "or a finite-difference scheme "
                    f"in {FINITE_DIFFERENCE_METHODS}, got {c_jac!r}."
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
            if callable(c_jac):
                # The Jacobian at x0 gives the sparsity pattern;
                # it is reused for the first evaluation by Uno.
                c_jac = _JacobianBlock(c_jac(x0.copy()), x0, m_i, i, jacobian=c_jac)
            else:
                # The finite differences give a dense Jacobian.
                c_fun = _LastValue(c_fun)
                c_jac = _JacobianBlock(
                    None,
                    x0,
                    m_i,
                    i,
                    jacobian=_create_finite_difference_jacobian(con.fun, c_fun, c_jac),
                )
            normalized.append(_Constraint(c_fun, c_jac, c_hess, c_lb, c_ub, m_i))
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
            normalized.append(_Constraint(c_fun, c_jac, c_hess, c_lb, c_ub, m_i))
        elif isinstance(con, dict):
            c_fun_raw = con["fun"]
            c_jac_raw = con.get("jac")
            if c_jac_raw is None or c_jac_raw is False:
                c_jac_raw = "2-point"
            if not (callable(c_jac_raw) or _is_finite_difference_method(c_jac_raw)):
                raise TypeError(
                    f"The 'jac' of the constraint {i} must be a callable "
                    "or a finite-difference scheme "
                    f"in {FINITE_DIFFERENCE_METHODS}, got {c_jac_raw!r}."
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
            if callable(c_jac_raw):
                c_jac = lambda x, _j=c_jac_raw, _a=c_args: _j(x, *_a)
                # The Jacobian at x0 gives the sparsity pattern;
                # it is reused for the first evaluation by Uno.
                c_jac = _JacobianBlock(c_jac(x0.copy()), x0, m_i, i, jacobian=c_jac)
            else:
                # The finite differences give a dense Jacobian.
                c_fun = _LastValue(c_fun)
                c_jac = _JacobianBlock(
                    None,
                    x0,
                    m_i,
                    i,
                    jacobian=_create_finite_difference_jacobian(
                        lambda x, _f=c_fun_raw, _a=c_args: _f(x, *_a), c_fun, c_jac_raw
                    ),
                )
            c_hess = None  # Dict constraints don't support Hessians
            normalized.append(_Constraint(c_fun, c_jac, c_hess, c_lb, c_ub, m_i))
        else:
            raise TypeError(f"Unsupported constraint type: {type(con)}")

        if callable(hess) and c_hess is None and not isinstance(
            con, LinearConstraint
        ):
            raise ValueError(
                f"The constraint {i} has no callable Hessian, "
                "which is required when 'hess' is provided; "
                "use a NonlinearConstraint with a callable 'hess'."
            )

    return normalized


class _MergedConstraints:
    """The normalized constraints merged into a single block of constraints for Uno."""

    def __init__(self, constraints: list[_Constraint]) -> None:
        """
        Args:
            constraints: The normalized constraints.
        """
        self.constraints = constraints
        """The normalized constraints."""

        self.size = sum(item.size for item in constraints)
        """The total number of constraints."""

        if self.size > 0:
            self.lower_bounds = np.concatenate(
                [item.lower_bound for item in constraints]
            )
            self.upper_bounds = np.concatenate(
                [item.upper_bound for item in constraints]
            )
        else:
            self.lower_bounds = np.array([])
            self.upper_bounds = np.array([])

        # Sparsity pattern (COO), concatenation of the patterns of the blocks
        if self.size > 0:
            self.nnz = sum(item.jacobian.nnz for item in constraints)
            row_offsets = np.cumsum([0] + [item.size for item in constraints[:-1]])
            self.row_indices = np.concatenate(
                [
                    item.jacobian.rows + offset
                    for item, offset in zip(constraints, row_offsets)
                ]
            ).tolist()
            self.col_indices = np.concatenate(
                [item.jacobian.cols for item in constraints]
            ).tolist()
        else:
            self.nnz = 0
            self.row_indices = []
            self.col_indices = []

    def evaluate(self, x: np.ndarray, constraint_values: np.ndarray) -> None:
        """Evaluate the constraints.

        Args:
            x: The point at which the constraints are evaluated.
            constraint_values: The output vector of the constraint values.
        """
        offset = 0
        for c_fun_i, _, _, _, _, m_i in self.constraints:
            constraint_values[offset:offset + m_i] = np.asarray(
                c_fun_i(x), dtype=float
            ).ravel()
            offset += m_i

    def evaluate_jacobian(self, x: np.ndarray, jacobian_values: np.ndarray) -> None:
        """Evaluate the Jacobian of the constraints.

        The values of each block are written in the order of the sparsity pattern,
        without building the dense merged Jacobian.

        Args:
            x: The point at which the Jacobian is evaluated.
            jacobian_values: The output vector of the Jacobian values.
        """
        start = 0
        for _, block, _, _, _, _ in self.constraints:
            block.write(x, jacobian_values[start:start + block.nnz])
            start += block.nnz
