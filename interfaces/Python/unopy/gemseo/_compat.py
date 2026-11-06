# Copyright (c) 2021-2026 IRT Saint Exupéry, https://www.irt-saintexupery.com
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""Compatibility layer between GEMSEO 6 and GEMSEO 7.

GEMSEO 7 moved most modules (``gemseo.algos.*`` to ``gemseo.optimization.*``
and ``gemseo.core.algorithm.*``) and changed a few library APIs.
The GEMSEO 7 names are tried first,
as GEMSEO 7 still provides the GEMSEO 6 ones as deprecated aliases.
"""

from __future__ import annotations

from typing import TYPE_CHECKING
from typing import Any

try:
    from gemseo.optimization.core.base_gradient_based_algorithm_settings import (
        BaseGradientBasedAlgorithmSettings,
    )
    from gemseo.optimization.core.base_optimization_library import (
        BaseOptimizationLibrary,
    )
    from gemseo.optimization.core.base_optimization_library import (
        OptimizationAlgorithmDescription,
    )
    from gemseo.optimization.core.base_optimizer_settings import BaseOptimizerSettings
    from gemseo.space.util import get_value_and_bounds
    from gemseo.util.constant import infinite_int as MAX_INT  # noqa: N812

    IS_GEMSEO_6 = False
except ImportError:  # GEMSEO 6
    from gemseo.algos.design_space_utils import get_value_and_bounds
    from gemseo.algos.opt.base_gradient_based_algorithm_settings import (
        BaseGradientBasedAlgorithmSettings,
    )
    from gemseo.algos.opt.base_optimization_library import BaseOptimizationLibrary
    from gemseo.algos.opt.base_optimization_library import (
        OptimizationAlgorithmDescription,
    )
    from gemseo.algos.opt.base_optimizer_settings import BaseOptimizerSettings
    from gemseo.utils.constants import C_LONG_MAX as MAX_INT  # noqa: N811

    IS_GEMSEO_6 = True

if TYPE_CHECKING:
    from gemseo.typing import RealArray

__all__ = [
    "IS_GEMSEO_6",
    "MAX_INT",
    "SETTINGS_FIELD",
    "BaseGradientBasedAlgorithmSettings",
    "BaseOptimizationLibrary",
    "BaseOptimizerSettings",
    "OptimizationAlgorithmDescription",
    "filter_settings",
    "get_initial_point_and_bounds",
]

SETTINGS_FIELD = "Settings" if IS_GEMSEO_6 else "settings_class"
"""The name of the field of an algorithm description holding its settings class."""


def get_initial_point_and_bounds(
    library: BaseOptimizationLibrary, problem: Any
) -> tuple[RealArray, RealArray, RealArray]:
    """Return the initial point and the bounds in the coordinates of the solver.

    Args:
        library: The optimization library solving the problem.
        problem: The optimization problem.

    Returns:
        The initial point, the lower bounds and the upper bounds.
    """
    if IS_GEMSEO_6:
        return get_value_and_bounds(
            problem.design_space, library._settings.normalize_design_space
        )
    # The problem works in the coordinates of the solver,
    # so its input space is already normalized if required.
    return get_value_and_bounds(problem.input_space)


def filter_settings(library: BaseOptimizationLibrary) -> dict[str, Any]:
    """Return the settings specific to the library, without the unset ones.

    Args:
        library: The optimization library.

    Returns:
        The settings that are neither GEMSEO optimizer settings nor ``None``.
    """
    if IS_GEMSEO_6:
        settings = library._filter_settings(
            library._settings.model_dump(), BaseOptimizerSettings
        )
    else:
        settings = library._filter_settings()
    return {name: value for name, value in settings.items() if value is not None}
