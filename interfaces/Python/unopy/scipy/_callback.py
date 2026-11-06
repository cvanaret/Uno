# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""The user callback called at each accepted iterate."""

from __future__ import annotations

import inspect
from typing import Any
from typing import Callable
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import OptimizeResult

if TYPE_CHECKING:
    from unopy.scipy._guard import _EvaluationErrorGuard
    from unopy.scipy._objective import _Objective

STOP_ITERATION_MESSAGE = "`callback` raised `StopIteration`."
"""The message of the result when the callback raised ``StopIteration``."""


def _uses_intermediate_result(callback: Callable[..., Any]) -> bool:
    """Return whether a callback takes an ``OptimizeResult``, as in SciPy.

    Args:
        callback: The callback.

    Returns:
        Whether the only parameter of the callback is named ``intermediate_result``.
    """
    try:
        parameters = inspect.signature(callback).parameters
    except (TypeError, ValueError):
        return False
    return set(parameters) == {"intermediate_result"}


class _IterationCallback:
    """The user callback called by Uno at each accepted iterate.

    The callback is called with either an ``OptimizeResult``
    (if its only parameter is named ``intermediate_result``) or a copy of ``x``.
    If it raises ``StopIteration``, Uno is stopped and the exception is not re-raised.
    Any other exception is handled as an exception raised by a user function:
    Uno is stopped and the exception is re-raised after the optimization.
    """

    def __init__(
        self,
        callback: Callable[..., Any],
        n: int,
        objective: _Objective,
        guard: _EvaluationErrorGuard,
    ) -> None:
        """
        Args:
            callback: The user callback.
            n: The number of variables.
            objective: The objective, giving its last value.
            guard: The guard of the user functions.
        """
        self._callback = callback
        self._uses_intermediate_result = _uses_intermediate_result(callback)
        self._n = n
        self._objective = objective
        self._guard = guard
        self.stopped = False
        """Whether the callback raised ``StopIteration``."""

        self.calls = 0
        """The number of calls to the callback."""

    def notify(
        self,
        primals: np.ndarray,
        lower_bound_multipliers: np.ndarray,
        upper_bound_multipliers: np.ndarray,
        constraint_multipliers: np.ndarray,
        objective_multiplier: float,
        primal_feasibility: float,
        stationarity: float,
        complementarity: float,
    ) -> None:
        """Call the user callback at an accepted iterate.

        Uno may pass the primals and the bound multipliers of a reformulated problem
        (e.g. with slack variables): only the first ``n`` entries are kept.
        Uno passes the multipliers with its internal sign convention,
        i.e. the opposite of the one of the result.

        Args:
            primals: The primal variables.
            lower_bound_multipliers: The multipliers of the lower bounds.
            upper_bound_multipliers: The multipliers of the upper bounds.
            constraint_multipliers: The multipliers of the constraints.
            objective_multiplier: The objective multiplier.
            primal_feasibility: The primal feasibility residual.
            stationarity: The stationarity residual.
            complementarity: The complementarity residual.
        """
        if self.stopped or self._guard.exception is not None:
            return
        try:
            x = np.array(primals[: self._n], dtype=float)
            if self._uses_intermediate_result:
                intermediate_result = OptimizeResult(
                    x=x,
                    fun=self._objective.value_at(x),
                    nit=self.calls + 1,
                    constraint_multipliers=-np.array(
                        constraint_multipliers, dtype=float
                    ),
                    lower_bound_multipliers=-np.array(
                        lower_bound_multipliers[: self._n], dtype=float
                    ),
                    upper_bound_multipliers=-np.array(
                        upper_bound_multipliers[: self._n], dtype=float
                    ),
                    objective_multiplier=float(objective_multiplier),
                    primal_feasibility=float(primal_feasibility),
                    stationarity=float(stationarity),
                    complementarity=float(complementarity),
                )
                self.calls += 1
                self._callback(intermediate_result=intermediate_result)
            else:
                self.calls += 1
                self._callback(x)
        except StopIteration:
            self.stopped = True
        except BaseException as exception:
            self._guard.record(exception)

    def terminate(self, *arguments: Any) -> bool:
        """Return whether Uno must stop.

        Args:
            *arguments: The arguments passed by Uno (iterate, multipliers and residuals).

        Returns:
            Whether the callback raised ``StopIteration``
            or a user function raised an exception.
        """
        return self.stopped or self._guard.terminate(*arguments)
