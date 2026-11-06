# Copyright (c) 2021-2026 IRT Saint Exupéry, https://www.irt-saintexupery.com
# Licensed under the MIT license. See LICENSE file in the project directory for details.
# Author: François Gallard
"""The library of Uno constrained gradient-based optimization algorithms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any
from typing import ClassVar
from typing import TYPE_CHECKING

from gemseo.algos.design_space_utils import get_value_and_bounds
from gemseo.algos.opt.base_optimization_library import BaseOptimizationLibrary
from gemseo.algos.opt.base_optimization_library import (
    OptimizationAlgorithmDescription,
)
from gemseo.algos.opt.base_optimizer_settings import BaseOptimizerSettings
from gemseo.utils.constants import C_LONG_MAX
from numpy import atleast_1d
from numpy import atleast_2d
from numpy import concatenate
from numpy import isfinite
from numpy import real
from numpy import vstack
from scipy.sparse import csr_array
from scipy.sparse import issparse
from scipy.sparse import vstack as sparse_vstack

from unopy.gemseo.settings.base_uno_settings import BaseUnoSettings
from unopy.gemseo.settings.filter_slp import UNO_Filter_SLP_Settings
from unopy.gemseo.settings.filter_sqp import UNO_Filter_SQP_Settings
from unopy.gemseo.settings.funnel_sqp import UNO_Funnel_SQP_Settings
from unopy.gemseo.settings.ipopt import UNO_IPOPT_Settings
from unopy.scipy import minimize

if TYPE_CHECKING:
    from collections.abc import Sequence

    from gemseo.algos.optimization_problem import OptimizationProblem
    from gemseo.core.mdo_functions.mdo_function import MDOFunction
    from gemseo.typing import NumberArray


def _concatenate_constraints(
    constraints: Sequence[MDOFunction], constraint_type: str
) -> dict[str, Any]:
    """Merge constraints of the same type into a single SciPy constraint.

    GEMSEO's ``Concatenate`` function is not used
    because it stacks the Jacobians with ``numpy.vstack``,
    which does not support sparse Jacobians,
    and because it would add a function wrapper.

    Args:
        constraints: The constraints of the same type.
        constraint_type: The type of the constraints, either "eq" or "ineq".

    Returns:
        The SciPy dictionary constraint concatenating the values of the constraints
        and stacking their Jacobians in the same order;
        the Jacobian is sparse if one of the Jacobians is sparse.
    """

    def evaluate(x: NumberArray) -> NumberArray:
        return concatenate([
            atleast_1d(constraint.evaluate(x)) for constraint in constraints
        ])

    def compute_jacobian(x: NumberArray) -> Any:
        jacobians = [constraint.jac(x) for constraint in constraints]
        if any(issparse(jacobian) for jacobian in jacobians):
            return sparse_vstack(
                [
                    jacobian if issparse(jacobian) else csr_array(atleast_2d(jacobian))
                    for jacobian in jacobians
                ],
                format="csr",
            )
        return vstack([atleast_2d(jacobian) for jacobian in jacobians])

    return {"type": constraint_type, "fun": evaluate, "jac": compute_jacobian}


def _merge_constraints(constraints: Sequence[MDOFunction]) -> list[dict[str, Any]]:
    """Merge the constraints into one equality and one inequality SciPy constraint.

    The order of the constraints of the same type is preserved.
    The functions of GEMSEO are evaluated with their own cache (the database),
    so a merged constraint does not evaluate a function twice at the same point.

    Args:
        constraints: The constraints with the sign expected by SciPy.

    Returns:
        At most two SciPy dictionary constraints,
        the equality one first.
    """
    merged_constraints = []
    for constraint_type in ("eq", "ineq"):
        constraints_of_type = [
            constraint for constraint in constraints
            if constraint.f_type == constraint_type
        ]
        if constraints_of_type:
            merged_constraints.append(
                _concatenate_constraints(constraints_of_type, constraint_type)
            )
    return merged_constraints


@dataclass
class UnoAlgorithmDescription(OptimizationAlgorithmDescription):
    """The description of the Uno constrained gradient-based optimization library."""

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

    Settings: type[BaseUnoSettings] = BaseUnoSettings
    """The option validation model for Uno optimization library."""

    website: str = "https://unosolver.readthedocs.io/en/latest/"
    """The website of the wrapped library or algorithm."""


class UnoOpt(BaseOptimizationLibrary[BaseUnoSettings]):
    """The library of Uno optimization algorithms."""

    ALGORITHM_INFOS: ClassVar[dict[str, UnoAlgorithmDescription]] = {
        "UNO_Filter_SQP": UnoAlgorithmDescription(
            algorithm_name="UNO_Filter_SQP",
            description=(
                "Sequential Quadratic Programming (SQP) "
                "implemented in the Uno library"
            ),
            internal_algorithm_name="filtersqp",
            Settings=UNO_Filter_SQP_Settings,
        ),
        "UNO_Filter_SLP": UnoAlgorithmDescription(
            algorithm_name="UNO_Filter_SLP",
            description=(
                "Sequential Linear Programming (SLP) implemented in the Uno library"
            ),
            internal_algorithm_name="filterslp",
            Settings=UNO_Filter_SLP_Settings,
        ),
        "UNO_Funnel_SQP": UnoAlgorithmDescription(
            algorithm_name="UNO_Funnel_SQP",
            description=(
                "Funnel Sequential Quadratic Programming (SQP) "
                "implemented in the Uno library"
            ),
            internal_algorithm_name="funnelsqp",
            Settings=UNO_Funnel_SQP_Settings,
        ),
        "UNO_IPOPT": UnoAlgorithmDescription(
            algorithm_name="UNO_IPOPT",
            description=(
                "Interior Point Optimization (IPOPT) implemented in the Uno library"
            ),
            internal_algorithm_name="ipopt",
            Settings=UNO_IPOPT_Settings,
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

        # Get the constraints in SciPy format,
        # merged into one equality and one inequality constraint.
        scipy_constraints = _merge_constraints(
            list(self._get_right_sign_constraints(problem))
        )

        # Filter settings to get only the Uno ones,
        # i.e. remove the GEMSEO ones and the Uno ones left to their default value.
        settings_ = self._filter_settings(
            self._settings.model_dump(exclude_none=True), BaseOptimizerSettings
        )

        # Deactivate stopping criteria which are handled by GEMSEO
        tolerance = 0.0
        # C_LONG_MAX is the largest 32-bit integer, i.e. the largest uno_int.
        settings_["max_iterations"] = C_LONG_MAX

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

        return str(opt_result.message), int(opt_result.status)
