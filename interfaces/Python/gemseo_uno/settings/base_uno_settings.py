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
"""Settings for the SciPy COBYQA algorithm."""

from __future__ import annotations

from gemseo.algos.opt.base_gradient_based_algorithm_settings import (
    BaseGradientBasedAlgorithmSettings,
)
from gemseo.algos.opt.base_optimizer_settings import BaseOptimizerSettings
from pydantic import PositiveFloat  # noqa:TC002


class UNO_Settings(
    BaseOptimizerSettings, BaseGradientBasedAlgorithmSettings
):  # noqa: N801
    """Settings for the UNO algorithm."""

    _TARGET_CLASS_NAME = "UNO"
    initial_tr_radius: PositiveFloat | None = None
    primal_tolerance: float | None = None
    dual_tolerance: float | None = None
    loose_primal_tolerance: float | None = None
    loose_dual_tolerance: float | None = None
    loose_tolerance_consecutive_iteration_threshold: int | None = None
    time_limit: float | None = None
    print_solution: bool | None = None
    unbounded_objective_threshold: float | None = None
    enforce_linear_constraints: bool | None = None
    logger: str | None = None
    constraint_relaxation_strategy: str | None = None
    inequality_handling_method: str | None = None
    globalization_mechanism: str | None = None
    globalization_strategy: str | None = None
    hessian_model: str | None = None
    inertia_correction_strategy: str | None = None
    scale_functions: bool | None = None
    function_scaling_threshold: float | None = None
    function_scaling_factor: float | None = None
    scale_residuals: bool | None = None
    progress_norm: str | None = None
    residual_norm: str | None = None
    residual_scaling_threshold: float | None = None
    protect_actual_reduction_against_roundoff: bool | None = None
    print_subproblem: bool | None = None
    armijo_decrease_fraction: float | None = None
    armijo_tolerance: float | None = None
    switching_delta: float | None = None
    switching_infeasibility_exponent: float | None = None
    filter_type: str | None = None
    filter_beta: float | None = None
    filter_gamma: float | None = None
    filter_ubd: float | None = None
    filter_fact: float | None = None
    filter_capacity: int | None = None
    filter_sufficient_infeasibility_decrease_factor: float | None = None
    nonmonotone_filter_number_dominated_entries: int | None = None
    funnel_kappa: float | None = None
    funnel_beta: float | None = None
    funnel_gamma: float | None = None
    funnel_ubd: float | None = None
    funnel_fact: float | None = None
    funnel_update_strategy: int | None = None
    funnel_require_acceptance_wrt_current_iterate: bool | None = None
    LS_backtracking_ratio: float | None = None
    LS_min_step_length: float | None = None
    LS_scale_duals_with_step_length: bool | None = None
    regularization_failure_threshold: float | None = None
    regularization_initial_value: float | None = None
    regularization_increase_factor: float | None = None
    primal_regularization_initial_factor: float | None = None
    dual_regularization_fraction: float | None = None
    primal_regularization_lb: float | None = None
    primal_regularization_decrease_factor: float | None = None
    primal_regularization_fast_increase_factor: float | None = None
    primal_regularization_slow_increase_factor: float | None = None
    threshold_unsuccessful_attempts: int | None = None
    TR_radius: float | None = None
    TR_increase_factor: float | None = None
    TR_decrease_factor: float | None = None
    TR_aggressive_decrease_factor: float | None = None
    TR_activity_tolerance: float | None = None
    TR_min_radius: float | None = None
    TR_radius_reset_threshold: float | None = None
    switch_to_optimality_requires_linearized_feasibility: bool | None = None
    l1_constraint_violation_coefficient: float | None = None
    barrier_initial_parameter: float | None = None
    barrier_default_multiplier: float | None = None
    barrier_tau_min: float | None = None
    barrier_k_sigma: float | None = None
    barrier_smax: float | None = None
    barrier_k_mu: float | None = None
    barrier_theta_mu: float | None = None
    barrier_k_epsilon: float | None = None
    barrier_update_fraction: float | None = None
    barrier_regularization_exponent: float | None = None
    barrier_small_direction_factor: float | None = None
    barrier_push_variable_to_interior_k1: float | None = None
    barrier_push_variable_to_interior_k2: float | None = None
    barrier_damping_factor: float | None = None
    least_square_multiplier_max_norm: float | None = None
    BQPD_kmax_heuristic: str | None = None
    QP_solver: str | None = None
    LP_solver: str | None = None
    linear_solver: str | None = None
