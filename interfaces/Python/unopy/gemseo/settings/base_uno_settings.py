# Copyright (c) 2021-2026 IRT Saint Exupéry, https://www.irt-saintexupery.com
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""The base settings of the Uno algorithms."""

from __future__ import annotations

from unopy.gemseo._compat import BaseGradientBasedAlgorithmSettings
from unopy.gemseo._compat import BaseOptimizerSettings


class BaseUnoSettings(BaseOptimizerSettings, BaseGradientBasedAlgorithmSettings):
    """The settings common to the Uno algorithms.

    The Uno options are optional fields named after the Uno options
    (see ``uno/options/Options.cpp``); only those set by the user are passed to Uno.
    The following Uno options are not exposed:

    - ``preset``, set from the algorithm name,
    - ``max_iterations``, set to the largest value
      so that GEMSEO stops the driver with ``max_iter``.

    The wrapper sets ``primal_tolerance`` and ``dual_tolerance`` to 0
    so that GEMSEO stops the driver with its own criteria;
    setting them here overrides these zero values.
    """

    primal_tolerance: float | None = None
    dual_tolerance: float | None = None
    loose_primal_tolerance: float | None = None
    loose_dual_tolerance: float | None = None
    loose_tolerance_iteration_threshold: int | None = None
    time_limit: float | None = None
    print_solution: bool | None = None
    unbounded_objective_threshold: float | None = None
    logger: str | None = None
    constraint_relaxation_strategy: str | None = None
    inequality_handling_method: str | None = None
    globalization_mechanism: str | None = None
    globalization_strategy: str | None = None
    hessian_model: str | None = None
    inertia_correction_strategy: str | None = None
    progress_norm: str | None = None
    residual_norm: str | None = None
    residual_scaling_threshold: float | None = None
    protect_actual_reduction_against_roundoff: bool | None = None
    protected_actual_reduction_macheps_coefficient: float | None = None
    print_subproblem: bool | None = None
    armijo_decrease_fraction: float | None = None
    armijo_tolerance: float | None = None
    switching_delta: float | None = None
    switching_infeasibility_exponent: float | None = None
    sufficient_infeasibility_decrease_ratio: float | None = None
    filter_type: str | None = None
    filter_beta: float | None = None
    filter_gamma: float | None = None
    filter_ubd: float | None = None
    filter_fact: float | None = None
    filter_capacity: int | None = None
    filter_sufficient_infeasibility_decrease_factor: float | None = None
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
    quasi_newton_memory_size: int | None = None
    quasi_newton_delta_lower_bound: float | None = None
    quasi_newton_delta_upper_bound: float | None = None
    LBFGS_max_skips_before_reset: int | None = None
    LSR1_pivot_max_magnitude: float | None = None
    regularization_failure_threshold: float | None = None
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
    barrier_function: str | None = None
    barrier_initial_parameter: float | None = None
    barrier_default_multiplier: float | None = None
    barrier_tau_min: float | None = None
    barrier_k_sigma: float | None = None
    barrier_k_mu: float | None = None
    barrier_theta_mu: float | None = None
    barrier_k_epsilon: float | None = None
    barrier_update_fraction: float | None = None
    barrier_regularization_exponent: float | None = None
    barrier_small_direction_factor: float | None = None
    barrier_push_variable_to_interior_k1: float | None = None
    barrier_push_variable_to_interior_k2: float | None = None
    barrier_damping_factor: float | None = None
    barrier_small_infeasibility_factor: float | None = None
    least_square_multiplier_max_norm: float | None = None
    BQPD_kmax_heuristic: str | None = None
    QP_solver: str | None = None
    LP_solver: str | None = None
    linear_solver: str | None = None
    option_file: str | None = None
    write_solution_to_file: bool | None = None
    SOC_max_iterations: int | None = None
    SOC_infeasibility_fraction: float | None = None
