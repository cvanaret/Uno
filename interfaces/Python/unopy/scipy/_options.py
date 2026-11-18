# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""The setting of the Uno options.

The available options and their types are listed in the documentation of ``minimize``.
"""

from __future__ import annotations

from typing import Any
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Mapping

    import unopy

_MAX_ITERATIONS_ALIASES = ("maxiter", "max_iter")
"""The aliases of the Uno option "max_iterations"."""


def _set_option(uno_solver: unopy.UnoSolver, name: str, value: Any) -> None:
    """Set an option of Uno.

    The bindings choose the type of the option from the Python type of the value,
    so an integer value of a floating-point option is converted to a float.

    Args:
        uno_solver: The Uno solver.
        name: The name of the option.
        value: The value of the option.

    Raises:
        TypeError: If the option does not accept the type of the value.
    """
    try:
        uno_solver.set_option(name, value)
    except TypeError as error:
        if not isinstance(value, (bool, np.bool_)) and isinstance(
            value, (int, np.integer)
        ):
            try:
                uno_solver.set_option(name, float(value))
                return
            except TypeError:
                pass
        # Uno's message names the type of the value, not the type of the option.
        raise TypeError(
            f"The Uno option {name!r} does not accept the value {value!r} "
            f"of type {type(value).__name__} ({error})."
        ) from None


def _set_options(
    uno_solver: unopy.UnoSolver, tol: float | None, options: Mapping[str, Any]
) -> None:
    """Set the tolerances and the options of Uno.

    Args:
        uno_solver: The Uno solver.
        tol: The tolerance used for both the primal and dual tolerances, if any.
        options: The Uno options,
            where "maxiter" and "max_iter" are aliases of "max_iterations".
    """
    if tol is not None:
        uno_solver.set_option("primal_tolerance", float(tol))
        uno_solver.set_option("dual_tolerance", float(tol))
    for key, value in options.items():
        if key in _MAX_ITERATIONS_ALIASES:
            key = "max_iterations"
        _set_option(uno_solver, key, value)
