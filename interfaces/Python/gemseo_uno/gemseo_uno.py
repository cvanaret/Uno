# Copyright 2021 IRT Saint Exupéry, https://www.irt-saintexupery.com
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License version 3 as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# Contributors:
#    INITIAL AUTHORS - initial API and implementation and/or initial
#                           documentation
#        :author: François Gallard
#    OTHER AUTHORS   - MACROSCOPIC CHANGES
"""The library of Uno constrained gradient-based optimization algorithms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any
from typing import ClassVar
from typing import TYPE_CHECKING

from gemseo.algos.design_space_utils import get_value_and_bounds
from gemseo.algos.opt.base_optimization_library import (
    OptimizationAlgorithmDescription,
    BaseOptimizationLibrary,
)
from numpy import isfinite
from numpy import real

from gemseo_uno.settings.base_uno_settings import UNO_Settings
from scipy_interface.scipy_uno import minimize

if TYPE_CHECKING:
    from gemseo.algos.optimization_problem import OptimizationProblem


@dataclass
class UnoAlgorithmDescription(OptimizationAlgorithmDescription):
    """The description of the Uno constrained grdient based optimization library."""

    library_name: str = "Uno"
    """The library name."""

    handle_equality_constraints: bool = True
    """Whether the optimization algorithm handles equality constraints."""

    handle_inequality_constraints: bool = True
    """Whether the optimization algorithm handles inequality constraints."""

    positive_constraints: bool = True
    """Whether the optimization algorithm requires positive constraints."""

    require_gradient: bool = True
    """Whether the optimization algorithm requires the gradient."""

    Settings: type[UNO_Settings] = UNO_Settings
    """The option validation model for Uno optimization library."""

    website: str = "https://unosolver.readthedocs.io/en/latest/"
    """The website of the wrapped library or algorithm."""


class UnoOpt(BaseOptimizationLibrary[UNO_Settings]):
    """The library of Uno optimization algorithms."""

    ALGORITHM_INFOS: ClassVar[dict[str, UnoAlgorithmDescription]] = {
        "UNO_Filter_SQP": UnoAlgorithmDescription(
            algorithm_name="UNO_Filter_SQP",
            description=(
                "Sequential Quadratic Programming (SQP) "
                "implemented in the Uno library"
            ),
            internal_algorithm_name="filtersqp",
            Settings=UNO_Settings,
        ),
        "UNO_Filter_SLP": UnoAlgorithmDescription(
            algorithm_name="UNO_Filter_SLP",
            description=(
                "Sequential Linear Programming (SLP) " "implemented in the Uno library"
            ),
            internal_algorithm_name="filterslp",
            Settings=UNO_Settings,
        ),
        "UNO_Funnel_SQP": UnoAlgorithmDescription(
            algorithm_name="UNO_Funnel_SQP",
            description=(
                "Funnel Sequential Quadratic Programming (SQP)"
                "implemented in the Uno library"
            ),
            internal_algorithm_name="funnelsqp",
            Settings=UNO_Settings,
        ),
        "UNO_IPOPT": UnoAlgorithmDescription(
            algorithm_name="UNO_IPOPT",
            description=(
                "Interior Point Optimization (IPOPT)" "implemented in the Uno library"
            ),
            internal_algorithm_name="ipopt",
            Settings=UNO_Settings,
        ),
    }

    def _run(self, problem: OptimizationProblem) -> tuple[str, Any]:
        # Get the normalized bounds:
        x_0, l_b, u_b = get_value_and_bounds(
            problem.design_space, self._settings.normalize_design_space
        )
        # Replace infinite values with None:
        l_b = [val if isfinite(val) else None for val in l_b]
        u_b = [val if isfinite(val) else None for val in u_b]
        bounds = list(zip(l_b, u_b, strict=False))

        # Get constraint in SciPy format
        scipy_constraints = [
            {
                "type": constraint.f_type,
                "fun": constraint.evaluate,
                "jac": constraint.jac,
            }
            for constraint in self._get_right_sign_constraints(problem)
        ]

        # Filter settings to get only the uno ones
        settings_ = self._filter_settings(self._settings.model_dump(), UNO_Settings)

        # Deactivate stopping criteria which are handled by GEMSEO
        tolerance = 0.0

        opt_result = minimize(
            fun=lambda x: real(problem.objective.evaluate(x)),
            jac=problem.objective.jac,
            x0=x_0,
            method=self.ALGORITHM_INFOS[self._algo_name].internal_algorithm_name,
            bounds=bounds,
            constraints=scipy_constraints,
            options=settings_,
            tol=tolerance,
        )

        return opt_result.message, opt_result.status
