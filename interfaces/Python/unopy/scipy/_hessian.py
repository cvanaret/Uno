# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""The Lagrangian Hessian operator."""

from __future__ import annotations

from typing import Any
from typing import Callable
from typing import TYPE_CHECKING

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator

if TYPE_CHECKING:
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
    """The Lagrangian Hessian operator computing Hessian-vector products.

    The Lagrangian Hessian
    ``obj_mult * hess(x) + sum_i constraint_hess_i(x, v_i)`` is never formed:
    the product with a vector is the sum of the products of its terms with this vector.
    Each term is either a dense array, a ``scipy.sparse`` matrix or array,
    or a ``scipy.sparse.linalg.LinearOperator``.
    The terms are computed once per point and reused
    by the successive Hessian-vector products computed by the subproblem solver:
    the objective Hessian is computed again when ``x`` changes
    and the Hessian of a constraint when ``x`` or its multipliers change.
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
        # Uno passes views on its own buffers: copy them to build the cache keys.
        x = np.array(x, dtype=float)
        multipliers = np.array(multipliers, dtype=float)
        vector = np.asarray(vector, dtype=float)

        # A term is cached only once computed, so that a failed evaluation is not cached.
        if not self._is_cached(self._objective_cache, x, None, evaluate_at_x):
            self.evaluations += 1
            term = _as_hessian_operator(self._hess(x, *self._args))
            self._objective_cache = (x, None, term)
        product = float(obj_mult) * np.asarray(
            self._objective_cache[2] @ vector, dtype=float
        ).ravel()

        offset = 0
        for index, constraint in enumerate(self._constraints):
            size = constraint.size
            if constraint.hessian is not None:
                multipliers_i = multipliers[offset:offset + size]
                cache = self._constraint_caches[index]
                if not self._is_cached(cache, x, multipliers_i, evaluate_at_x):
                    term = _as_hessian_operator(constraint.hessian(x, multipliers_i))
                    cache = self._constraint_caches[index] = (x, multipliers_i, term)
                product += np.asarray(cache[2] @ vector, dtype=float).ravel()
            offset += size

        result[:] = product
