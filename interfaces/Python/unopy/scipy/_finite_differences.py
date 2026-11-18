# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""The finite-difference approximation of first-order derivatives.

SciPy has no public function approximating the Jacobian of a vector-valued function
with central differences or complex steps
(``scipy.optimize.approx_fprime`` uses forward differences only),
so a small implementation is provided here.
The steps are the ones of SciPy: ``h = rel_step * sign(x) * max(1, |x|)``
with ``rel_step`` equal to ``sqrt(eps)`` for "2-point" and "cs"
and to ``eps ** (1/3)`` for "3-point", where ``sign(0) = 1``.
The bounds of the variables are not taken into account:
the function may be evaluated slightly outside them.
"""

from __future__ import annotations

from typing import Any
from typing import Callable

import numpy as np

FINITE_DIFFERENCE_METHODS = ("2-point", "3-point", "cs")
"""The finite-difference schemes."""

_EPS = np.finfo(float).eps

_RELATIVE_STEPS = {"2-point": _EPS**0.5, "3-point": _EPS ** (1 / 3), "cs": _EPS**0.5}
"""The relative steps of the finite-difference schemes."""


def _is_finite_difference_method(value: Any) -> bool:
    """Return whether a value is the name of a finite-difference scheme.

    Args:
        value: The value.

    Returns:
        Whether the value is the name of a finite-difference scheme.
    """
    return isinstance(value, str) and value in FINITE_DIFFERENCE_METHODS


def _approximate_jacobian(
    function: Callable[[np.ndarray], Any],
    x: np.ndarray,
    method: str,
    value: np.ndarray | None = None,
) -> np.ndarray:
    """Approximate the Jacobian of a function with finite differences.

    Args:
        function: The function ``function(x) -> array_like`` of size ``m``.
            With "cs", it must accept a complex ``x``.
        x: The point at which the Jacobian is approximated.
        method: The finite-difference scheme, either "2-point", "3-point" or "cs".
        value: The value of the function at ``x``, used by "2-point".
            If ``None``, the function is evaluated at ``x``.

    Returns:
        The approximated Jacobian of shape ``(m, n)``.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    sign = np.where(x >= 0, 1.0, -1.0)
    steps = _RELATIVE_STEPS[method] * sign * np.maximum(1.0, np.abs(x))
    # Use the steps exactly representable as the difference of two floats.
    steps = (x + steps) - x
    columns = []
    if method == "2-point":
        if value is None:
            value = function(x.copy())
        value = np.asarray(value, dtype=float).ravel()
        for j in range(n):
            x_j = x.copy()
            x_j[j] += steps[j]
            columns.append(
                (np.asarray(function(x_j), dtype=float).ravel() - value) / steps[j]
            )
    elif method == "3-point":
        for j in range(n):
            x_forward = x.copy()
            x_forward[j] += steps[j]
            x_backward = x.copy()
            x_backward[j] -= steps[j]
            columns.append(
                (
                    np.asarray(function(x_forward), dtype=float).ravel()
                    - np.asarray(function(x_backward), dtype=float).ravel()
                )
                / (2 * steps[j])
            )
    else:
        for j in range(n):
            x_j = x.astype(complex)
            x_j[j] += 1j * steps[j]
            columns.append(np.asarray(function(x_j), dtype=complex).ravel().imag / steps[j])
    return np.column_stack(columns) if columns else np.empty((0, 0))
