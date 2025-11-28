# Options

This document lists the configuration options that Uno accepts, their types and default values. Options are grouped by functional area. Defaults are taken from `uno/options/DefaultOptions.cpp`.

## Ingredients
- `constraint_relaxation_strategy` (string) - {feasibility_restoration}
- `inequality_handling_method` (string) - {inequality_constrained, interior_point}
- `hessian_model` (string) - {exact, identity, zero}
- `inertia_correction_strategy` (string) - {primal, primal_dual, none}
- `globalization_mechanism` (string) - {TR, LS}
- `globalization_strategy` (string) - {l1_merit, fletcher_filter_method, waechter_filter_method, funnel_method}

## Subproblem solver
- `QP_solver` (string) - if not provided, chosen automatically from available QP solvers if any
- `LP_solver` (string) - if not provided, chosen automatically from available LP solvers if any
- `linear_solver` (string) - if not provided, chosen automatically from available linear solvers if any

## Termination
- `primal_tolerance` (double) - 1e-8
- `dual_tolerance` (double) - 1e-8
- `loose_primal_tolerance` (double) - 1e-6
- `loose_dual_tolerance` (double) - 1e-6
- `loose_tolerance_consecutive_iteration_threshold` (integer) - 15
- `max_iterations` (integer) - 2000
- `time_limit` (double) - infinity
- `print_solution` (bool) - false
- `unbounded_objective_threshold` (double) - -1e20
- `enforce_linear_constraints` (bool) - false

## Statistics table ordering
- `statistics_major_column_order` (integer) - 1
- `statistics_minor_column_order` (integer) - 2
- `statistics_penalty_parameter_column_order` (integer) - 5
- `statistics_barrier_parameter_column_order` (integer) - 8
- `statistics_SOC_column_order` (integer) - 9
- `statistics_TR_radius_column_order` (integer) - 10
- `statistics_LS_step_length_column_order` (integer) - 10
- `statistics_restoration_phase_column_order` (integer) - 20
- `statistics_primal_regularization_column_order` (integer) - 21
- `statistics_funnel_width_column_order` (integer) - 25
- `statistics_step_norm_column_order` (integer) - 31
- `statistics_objective_column_order` (integer) - 100
- `statistics_primal_feasibility_column_order` (integer) - 101
- `statistics_dual_feasibility_column_order` (integer) - 102
- `statistics_stationarity_column_order` (integer) - 104
- `statistics_complementarity_column_order` (integer) - 105
- `statistics_status_column_order` (integer) - 200

## Main options
- `logger` (string) - {SILENT, DISCRETE, WARNING, INFO, DEBUG, DEBUG2, DEBUG3} (default "INFO")
- `scale_functions` (bool) - false
- `function_scaling_threshold` (double) - 100.0
- `function_scaling_factor` (double) - 100.0
- `scale_residuals` (bool) - true
- `progress_norm` (string) - {L1, L2, INF} (default "L1")
- `residual_norm` (string) - {L1, L2, INF} (default "INF")
- `residual_scaling_threshold` (double) - 100.0
- `protect_actual_reduction_against_roundoff` (bool) - false
- `print_subproblem` (bool) - false

## Globalization strategy options
- `armijo_decrease_fraction` (double) - 1e-4
- `armijo_tolerance` (double) - 1e-9

## Switching method options
- `switching_delta` (double) - 0.999
- `switching_infeasibility_exponent` (double) - 2

## Filter method options
- `filter_type` (string) - {standard, nonmonotone} (default "standard")
- `filter_beta` (double) - 0.999
- `filter_gamma` (double) - 0.001
- `filter_ubd` (double) - 1e2
- `filter_fact` (double) - 1.25
- `filter_capacity` (integer) - 50
- `filter_sufficient_infeasibility_decrease_factor` (double) - 0.9
- `nonmonotone_filter_number_dominated_entries` (integer) - 3

## Funnel options
- `funnel_kappa` (double) - 0.5
- `funnel_beta` (double) - 0.9999
- `funnel_gamma` (double) - 0.001
- `funnel_ubd` (double) - 1.0
- `funnel_fact` (double) - 1.5
- `funnel_update_strategy` (integer) - 1
- `funnel_require_acceptance_wrt_current_iterate` (bool) - false

## Line search options
- `LS_backtracking_ratio` (double) - 0.5
- `LS_min_step_length` (double) - 1e-12
- `LS_scale_duals_with_step_length` (bool) - true

## Regularization options
- `regularization_failure_threshold` (double) - 1e40
- `regularization_initial_value` (double) - 1e-4
- `regularization_increase_factor` (double) - 2
- `primal_regularization_initial_factor` (double) - 1e-4
- `dual_regularization_fraction` (double) - 1e-8
- `primal_regularization_lb` (double) - 1e-20
- `primal_regularization_decrease_factor` (double) - 3.0
- `primal_regularization_fast_increase_factor` (double) - 100.0
- `primal_regularization_slow_increase_factor` (double) - 8.0
- `threshold_unsuccessful_attempts` (integer) - 8

## Trust region options
- `TR_radius` (double) - 10.0
- `TR_increase_factor` (double) - 2
- `TR_decrease_factor` (double) - 2
- `TR_aggressive_decrease_factor` (double) - 4
- `TR_activity_tolerance` (double) - 1e-6
- `TR_min_radius` (double) - 1e-7
- `TR_radius_reset_threshold` (double) - 1e-4

## Feasibility restoration options
- `switch_to_optimality_requires_linearized_feasibility` (bool) - true
- `l1_constraint_violation_coefficient` (double) - 1

## Barrier subproblem options
- `barrier_function` (string) - "log"
- `barrier_initial_parameter` (double) - 0.1
- `barrier_default_multiplier` (double) - 1
- `barrier_tau_min` (double) - 0.99
- `barrier_k_sigma` (double) - 1e10
- `barrier_smax` (double) - 100
- `barrier_k_mu` (double) - 0.2
- `barrier_theta_mu` (double) - 1.5
- `barrier_k_epsilon` (double) - 10
- `barrier_update_fraction` (double) - 10
- `barrier_regularization_exponent` (double) - 0.25
- `barrier_small_direction_factor` (double) - 10.0
- `barrier_push_variable_to_interior_k1` (double) - 1e-2
- `barrier_push_variable_to_interior_k2` (double) - 1e-2
- `barrier_damping_factor` (double) - 1e-5
- `least_square_multiplier_max_norm` (double) - 1e3

## BQPD options
- `BQPD_kmax_heuristic` (string) - {"filtersqp", "minotaur"} (default "filtersqp")
