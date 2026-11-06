# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""The processing of the bounds of the variables."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import Bounds

if TYPE_CHECKING:
    from collections.abc import Iterable

Inf = float("inf")


def _process_bounds(
    bounds: Bounds | Iterable[tuple[float | None, float | None]], n: int
) -> tuple[np.ndarray, np.ndarray]:
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


def _get_variable_bounds(
    bounds: Bounds | Iterable[tuple[float | None, float | None]] | None, n: int
) -> tuple[np.ndarray | list[float], np.ndarray | list[float]]:
    """Return the lower and upper bounds of the variables.

    Args:
        bounds: Either a "Bounds" instance, a sequence of "(min, max)" pairs
            with "None" meaning no bound, or "None" for unbounded variables.
            It is not modified.
        n: The number of variables.

    Returns:
        The lower and upper bounds, each of size n.

    Raises:
        ValueError: If the bounds are ill-shaped or a lower bound is greater
            than the corresponding upper bound.
    """
    if bounds is not None:
        return _process_bounds(bounds, n)
    return [-Inf] * n, [Inf] * n
