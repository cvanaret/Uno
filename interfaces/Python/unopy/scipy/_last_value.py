# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""The storage of the last value of a function with its point."""

from __future__ import annotations

from typing import Any

import numpy as np


class _LastPointValue:
    """The last value of a function computed for Uno, stored with its point.

    The stored value is reused, e.g. by the finite differences ("2-point"),
    when the function is needed again at the same point.
    """

    def __init__(self) -> None:
        self._x = None
        """The last point at which the function has been computed for Uno."""

        self._value = None
        """The value of the function at ``_x``."""

    def _is_last_point(self, x: np.ndarray) -> bool:
        """Return whether a point is the last one at which the function was computed.

        Args:
            x: The point.

        Returns:
            Whether the point is the last one.
        """
        return self._x is not None and np.array_equal(x, self._x)

    def _store(self, x: np.ndarray, value: Any) -> None:
        """Store a value of the function with its point.

        Args:
            x: The point.
            value: The value of the function at ``x``.
        """
        # Uno passes a view on its own buffer: copy it.
        self._x = np.array(x, dtype=float)
        self._value = value
