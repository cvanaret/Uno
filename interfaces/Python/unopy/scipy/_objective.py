# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""The objective function and its gradient."""

from __future__ import annotations

from typing import Any
from typing import Callable

import numpy as np

from unopy.scipy._finite_differences import FINITE_DIFFERENCE_METHODS
from unopy.scipy._finite_differences import _approximate_jacobian
from unopy.scipy._finite_differences import _is_finite_difference_method
from unopy.scipy._last_value import _LastPointValue


class _Objective(_LastPointValue):
    """The objective function and its gradient computed for Uno.

    The last value of the objective computed for Uno is stored with its point,
    so that it is reused by the finite differences ("2-point"),
    by the intermediate results passed to the callback
    and, when ``jac=True``, by the gradient.
    """

    def __init__(
        self,
        fun: Callable[..., Any],
        jac: Callable[..., Any] | str | bool | None,
        args: tuple[Any, ...],
    ) -> None:
        """
        Args:
            fun: The objective function ``fun(x, *args) -> float``,
                or ``fun(x, *args) -> (float, array_like)`` if ``jac`` is ``True``.
            jac: The gradient of the objective function
                ``jac(x, *args) -> array_like``,
                ``True`` if ``fun`` returns the objective and its gradient,
                or a finite-difference scheme in ("2-point", "3-point", "cs");
                ``None`` and ``False`` mean "2-point".

        Raises:
            ValueError: If ``jac`` is not valid.
        """
        if jac is None or jac is False:
            jac = "2-point"
        if not (callable(jac) or jac is True or _is_finite_difference_method(jac)):
            raise ValueError(
                "The objective gradient 'jac' must be a callable, True, None "
                f"or a finite-difference scheme in {FINITE_DIFFERENCE_METHODS}, "
                f"got {jac!r}."
            )
        super().__init__()
        self._fun = fun
        self._jac = jac
        self._args = args
        self._gradient = None
        """The gradient of the objective at ``_x`` when ``jac`` is ``True``."""

        self.evaluations = 0
        """The number of calls to ``fun``, including the finite-difference ones."""

    def _call_fun(self, x: np.ndarray) -> Any:
        """Call the objective function and count the call.

        Args:
            x: The point.

        Returns:
            The output of the objective function.
        """
        self.evaluations += 1
        return self._fun(x, *self._args)

    @staticmethod
    def _to_scalar(value: Any, dtype: type = float) -> float | complex:
        """Convert the value of the objective to a scalar.

        Args:
            value: The value of the objective.
            dtype: The type of the scalar.

        Returns:
            The scalar value.

        Raises:
            ValueError: If the value is not of size 1.
        """
        value = np.asarray(value, dtype=dtype)
        if value.size != 1:
            raise ValueError(
                "The objective function must return a scalar, "
                f"got an array of shape {value.shape}."
            )
        return value.item()

    def evaluate(self, x: np.ndarray) -> float:
        """Evaluate the objective function for Uno.

        Args:
            x: The point.

        Returns:
            The value of the objective.
        """
        if self._jac is True and self._is_last_point(x):
            return self._value
        output = self._call_fun(x)
        if self._jac is True:
            value, gradient = output
            gradient = np.array(gradient, dtype=float).ravel()
        else:
            value, gradient = output, None
        value = self._to_scalar(value)
        self._store(x, value)
        self._gradient = gradient
        return value

    def value_at(self, x: np.ndarray) -> float:
        """Return the value of the objective, reusing the last one if possible.

        Args:
            x: The point.

        Returns:
            The value of the objective.
        """
        if self._is_last_point(x):
            return self._value
        return self.evaluate(x)

    def compute_gradient(self, x: np.ndarray, gradient: np.ndarray) -> None:
        """Compute the gradient of the objective function for Uno.

        Args:
            x: The point.
            gradient: The output vector of the gradient.
        """
        if callable(self._jac):
            gradient[:] = np.asarray(self._jac(x, *self._args), dtype=float).ravel()
        elif self._jac is True:
            if not self._is_last_point(x):
                self.evaluate(x)
            gradient[:] = self._gradient
        else:
            value = self.value_at(x) if self._jac == "2-point" else None
            gradient[:] = _approximate_jacobian(
                lambda x_: self._to_scalar(
                    self._call_fun(x_), complex if self._jac == "cs" else float
                ),
                x,
                self._jac,
                None if value is None else np.array([value]),
            ).ravel()
