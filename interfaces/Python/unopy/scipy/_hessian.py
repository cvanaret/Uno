# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""The Lagrangian Hessian, as an operator and as a sparse matrix."""

from __future__ import annotations

from typing import Any
from typing import Callable
from typing import TYPE_CHECKING

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator

if TYPE_CHECKING:
    from collections.abc import Iterator

    from unopy.scipy._constraints import _Constraint


def _as_hessian_operator(value: Any) -> np.ndarray | sp.spmatrix | LinearOperator:
    """Return a Hessian as an object supporting the product with a vector.

    A ``scipy.sparse`` matrix or array and a ``scipy.sparse.linalg.LinearOperator``
    are returned as is, so that no dense matrix is built;
    any other value is converted to a dense array of floats.

    Args:
        value: The Hessian.

    Returns:
        The Hessian supporting ``hessian @ vector``.
    """
    if sp.issparse(value) or isinstance(value, LinearOperator):
        return value
    return np.atleast_2d(np.asarray(value, dtype=float))


class _LagrangianHessianOperator:
    """The Lagrangian Hessian, as an operator and as a sparse matrix.

    The Lagrangian Hessian is ``obj_mult * hess(x) + sum_i constraint_hess_i(x, v_i)``.
    Each term is either a dense array, a ``scipy.sparse`` matrix or array,
    or a ``scipy.sparse.linalg.LinearOperator``.

    As an operator, it computes Hessian-vector products without forming the Hessian:
    the product with a vector is the sum of the products of its terms with this vector.
    This is used by the QP solver BQPD.

    As a sparse matrix, it computes the values of the lower triangle of the Hessian
    in a sparsity pattern fixed at ``x0`` (see ``create_sparsity_pattern``).
    This requires every term to be a dense array or a sparse matrix
    and is used by the direct linear solvers, e.g. with the IPOPT preset.

    The terms are computed once per point and shared by both representations:
    the objective Hessian is computed again when ``x`` changes
    and the Hessian of a constraint when ``x`` or its multipliers change.
    A term whose multipliers are all zero is not computed.
    """

    def __init__(
        self,
        hess: Callable[..., Any],
        args: tuple[Any, ...],
        constraints: list[_Constraint],
    ) -> None:
        """
        Args:
            hess: The Hessian of the objective function ``hess(x, *args)``.
            args: The extra arguments passed to ``hess``.
            constraints: The normalized constraints.
        """
        self._hess = hess
        self._args = args
        self._constraints = constraints
        # The last term of each Hessian and the point at which it was evaluated:
        # (x, None) for the objective, (x, multipliers) for a constraint.
        self._objective_cache = None
        self._constraint_caches = [None] * len(constraints)
        self.evaluations = 0
        """The number of evaluations of the user Hessian of the objective.

        Uno does not count the evaluations of a Hessian operator.
        """

    @staticmethod
    def _is_cached(
        cache: tuple[np.ndarray, np.ndarray | None, Any] | None,
        x: np.ndarray,
        multipliers: np.ndarray | None,
        evaluate_at_x: bool,
    ) -> bool:
        """Return whether a cached term can be reused.

        Args:
            cache: The cached ``(x, multipliers, term)``, if any.
            x: The point.
            multipliers: The multipliers of the term, if any.
            evaluate_at_x: Whether the point may have changed since the last call;
                if ``False``, the point is not compared with the cached one,
                but the multipliers still are.

        Returns:
            Whether the cached term can be reused.
        """
        if cache is None:
            return False
        if evaluate_at_x and not np.array_equal(x, cache[0]):
            return False
        return multipliers is None or np.array_equal(multipliers, cache[1])

    def _get_objective_term(self, x: np.ndarray, evaluate_at_x: bool) -> Any:
        """Return the Hessian of the objective, computed once per point.

        A term is cached only once computed, so that a failed evaluation is not cached.

        Args:
            x: The point, owned by the caller.
            evaluate_at_x: Whether the point may have changed since the last call.

        Returns:
            The Hessian of the objective.
        """
        if not self._is_cached(self._objective_cache, x, None, evaluate_at_x):
            self.evaluations += 1
            term = _as_hessian_operator(self._hess(x, *self._args))
            self._objective_cache = (x, None, term)
        return self._objective_cache[2]

    def _get_constraint_term(
        self,
        index: int,
        x: np.ndarray,
        multipliers: np.ndarray,
        evaluate_at_x: bool,
    ) -> Any:
        """Return the Hessian of a constraint, computed once per point and multipliers.

        Args:
            index: The index of the constraint.
            x: The point, owned by the caller.
            multipliers: The multipliers of the constraint, owned by the caller.
            evaluate_at_x: Whether the point may have changed since the last call.

        Returns:
            The Hessian of the constraint weighted by its multipliers.
        """
        cache = self._constraint_caches[index]
        if not self._is_cached(cache, x, multipliers, evaluate_at_x):
            hessian = self._constraints[index].hessian
            term = _as_hessian_operator(hessian(x, multipliers))
            cache = self._constraint_caches[index] = (x, multipliers, term)
        return cache[2]

    def _iterate_terms(
        self,
        x: np.ndarray,
        evaluate_at_x: bool,
        obj_mult: float,
        multipliers: np.ndarray,
    ) -> Iterator[tuple[float, Any]]:
        """Iterate over the terms of the Lagrangian Hessian with their weights.

        A term with zero multipliers is not computed,
        e.g. the objective term in the feasibility restoration phase.

        Args:
            x: The point.
            evaluate_at_x: Whether the point may have changed since the last call.
            obj_mult: The objective multiplier.
            multipliers: The constraint multipliers.

        Yields:
            The weight of a term and the term;
            the weight of a constraint term is 1
            since the constraint Hessian is already weighted by its multipliers.
        """
        # Uno passes views on its own buffers: copy them to build the cache keys.
        x = np.array(x, dtype=float)
        multipliers = np.array(multipliers, dtype=float)
        obj_mult = float(obj_mult)
        if obj_mult != 0.0:
            yield obj_mult, self._get_objective_term(x, evaluate_at_x)

        offset = 0
        for index, constraint in enumerate(self._constraints):
            size = constraint.size
            multipliers_i = multipliers[offset:offset + size]
            if constraint.hessian is not None and np.any(multipliers_i):
                yield 1.0, self._get_constraint_term(
                    index, x, multipliers_i, evaluate_at_x
                )
            offset += size

    def __call__(
        self,
        x: np.ndarray,
        evaluate_at_x: bool,
        obj_mult: float,
        multipliers: np.ndarray,
        vector: np.ndarray,
        result: np.ndarray,
    ) -> None:
        """Compute the product of the Lagrangian Hessian with a vector.

        Args:
            x: The point.
            evaluate_at_x: Whether the point may have changed since the last call.
            obj_mult: The objective multiplier.
            multipliers: The constraint multipliers.
            vector: The vector.
            result: The output vector of the Hessian-vector product.
        """
        vector = np.asarray(vector, dtype=float)
        product = np.zeros(len(vector))
        for weight, term in self._iterate_terms(
            x, evaluate_at_x, obj_mult, multipliers
        ):
            product += weight * np.asarray(term @ vector, dtype=float).ravel()
        result[:] = product

    def create_sparsity_pattern(self, x0: np.ndarray) -> bool:
        """Create the sparsity pattern of the Lagrangian Hessian as a matrix.

        The sparsity pattern is the union of the lower triangles of the terms at ``x0``:
        all the entries of a dense term and the stored entries of a sparse term.
        The Hessian of a constraint is evaluated with positive multipliers
        drawn at random, so that the Hessians of its components do not cancel out.
        The sparsity pattern is fixed for the whole optimization.

        Args:
            x0: The initial point.

        Returns:
            Whether the Lagrangian Hessian can be computed as a matrix,
            i.e. whether no term is a ``LinearOperator``.
        """
        x0 = np.array(x0, dtype=float)
        terms = [self._get_objective_term(x0, True)]
        generator = np.random.default_rng(0)
        for index, constraint in enumerate(self._constraints):
            if constraint.hessian is not None:
                multipliers = generator.uniform(1.0, 2.0, constraint.size)
                terms.append(self._get_constraint_term(index, x0, multipliers, True))
        if any(isinstance(term, LinearOperator) for term in terms):
            return False

        n = len(x0)
        keys = [
            self._get_lower_triangle(term, n, keep_zeros=True)[0] for term in terms
        ]
        self._hessian_keys = np.unique(np.concatenate(keys))
        self.hessian_rows, self.hessian_columns = np.divmod(self._hessian_keys, n)
        return True

    @staticmethod
    def _get_lower_triangle(
        term: np.ndarray | sp.spmatrix, n: int, keep_zeros: bool = False
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return the entries of the lower triangle of a term.

        Args:
            term: The term, either a dense array or a ``scipy.sparse`` matrix.
            n: The number of variables.
            keep_zeros: Whether to return all the entries of a dense term,
                including its zeros; the stored entries of a sparse term
                are always returned, including its explicit zeros.

        Returns:
            The linear indices ``row * n + column`` of the entries
            and their values.

        Raises:
            ValueError: If the term does not have the shape ``(n, n)``.
        """
        if term.shape != (n, n):
            raise ValueError(
                f"A Hessian must have the shape {(n, n)}, got {term.shape}."
            )
        if sp.issparse(term):
            term = sp.coo_array(term)
            rows = term.row.astype(np.int64)
            columns = term.col.astype(np.int64)
            data = term.data
        elif keep_zeros:
            rows, columns = np.tril_indices(n)
            data = term[rows, columns]
        else:
            rows, columns = np.nonzero(np.tril(term))
            data = term[rows, columns]
        in_lower_triangle = rows >= columns
        return (
            rows[in_lower_triangle].astype(np.int64) * n + columns[in_lower_triangle],
            np.asarray(data, dtype=float)[in_lower_triangle],
        )

    def evaluate_matrix(
        self,
        x: np.ndarray,
        obj_mult: float,
        multipliers: np.ndarray,
        hessian_values: np.ndarray,
    ) -> None:
        """Compute the values of the Lagrangian Hessian in its sparsity pattern.

        Args:
            x: The point.
            obj_mult: The objective multiplier.
            multipliers: The constraint multipliers.
            hessian_values: The output vector of the values of the lower triangle,
                in the order of the sparsity pattern.

        Raises:
            ValueError: If a term has a nonzero entry outside the sparsity pattern.
        """
        n = len(x)
        values = np.zeros(len(self._hessian_keys))
        for weight, term in self._iterate_terms(x, True, obj_mult, multipliers):
            keys, data = self._get_lower_triangle(term, n)
            positions = np.searchsorted(self._hessian_keys, keys)
            in_pattern = positions < len(self._hessian_keys)
            in_pattern[in_pattern] = (
                self._hessian_keys[positions[in_pattern]] == keys[in_pattern]
            )
            outside = ~in_pattern & (data != 0.0)
            if np.any(outside):
                row, column = divmod(int(keys[np.argmax(outside)]), n)
                raise ValueError(
                    f"A Hessian has a nonzero entry at ({row}, {column}) "
                    "outside the sparsity pattern of the Lagrangian Hessian, "
                    "which is the one of the Hessians at x0."
                )
            np.add.at(values, positions[in_pattern], weight * data[in_pattern])
        hessian_values[:] = values
