# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""The guard catching the exceptions raised by the user functions."""

from __future__ import annotations

from typing import Any
from typing import Callable


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

    def __init__(self) -> None:
        self.exception: BaseException | None = None
        """The first exception raised by a user function, if any."""

    def wrap_objective(
        self, callback: Callable[[Any], float]
    ) -> Callable[[Any], float]:
        """Wrap a callback returning the objective value.

        Args:
            callback: The callback ``callback(x) -> float``.

        Returns:
            The callback returning 0.0 instead of raising an exception.
        """

        def guarded_callback(x: Any) -> float:
            if self.exception is None:
                try:
                    return callback(x)
                except BaseException as exception:
                    self.exception = exception
            return 0.0

        return guarded_callback

    def wrap_inplace(
        self, callback: Callable[..., None], fill_with_zeros: bool = True
    ) -> Callable[..., None]:
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

        def guarded_callback(*arguments: Any) -> None:
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

    def record(self, exception: BaseException) -> None:
        """Store an exception raised by a user function if none has been stored.

        Args:
            exception: The exception.
        """
        if self.exception is None:
            self.exception = exception

    def terminate(self, *_: Any) -> bool:
        """Return whether Uno must stop because a user function raised an exception.

        Args:
            *_: The arguments passed by Uno (iterate, multipliers and residuals).

        Returns:
            Whether a user function raised an exception.
        """
        return self.exception is not None

    def raise_exception(self) -> None:
        """Re-raise the exception raised by a user function, if any.

        Raises:
            BaseException: The exception raised by a user function.
        """
        if self.exception is not None:
            raise self.exception
