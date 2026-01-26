# Options

This document lists the configuration options that Uno accepts, their types and default values. Options are grouped by functional area. Defaults are taken from `uno/options/DefaultOptions.cpp`.

## Ingredients

| Ingredient | Available strategies |
| :---       | :---                 |
| `constraint_relaxation_strategy` | `feasibility_restoration` |
| `inequality_handling_method` | `inequality_constrained`, `interior_point` |
| `hessian_model` | `exact`, `identity`, `zero` |
| `inertia_correction_strategy` | `primal`, `primal_dual`, `none` |
| `globalization_mechanism` | `TR`, `LS` |
| `globalization_strategy` | `l1_merit`, `fletcher_filter_method`, `waechter_filter_method`, `funnel_method` |

## Subproblem solver

| Ingredient | Available strategies |
| :---       | :---                 |
| `QP_solver` | `BQPD` (if available) |
| `LP_solver` | `BQPD`, `HiGHS` (if available) |
| `linear_solver` | `MA57`, `MA27`, `MUMPS` (if available) |

If not provided, the solver is chosen automatically from the available solvers (if any).

## Termination

| Option             | Type | Default | Description |
| :---               | :--- | :---    | :---        |
| `primal_tolerance` | double | 1e-8 | Tolerance on constraint violation |
| `dual_tolerance` | double | 1e-8 | Tolerance on stationarity and complementarity |
| `loose_primal_tolerance` | double | 1e-6 | Loose tolerance on constraint violation |
| `loose_dual_tolerance` | double | 1e-6 | Loose tolerance on stationarity and complementarity |
| `loose_tolerance_iteration_threshold` | integer | 15 | Number of iterations for the loose tolerance to apply |
| `max_iterations` | integer | 2000 | Maximum number of outer iterations |
| `time_limit` | double | infinity | Time limit |
| `print_solution` | bool | false | Whether the primal-dual solution is printed |
| `unbounded_objective_threshold` | double | -1e20 | Objective threshold under which the problem is declared unbounded |

## Main options

### String options
| Option             | Possible values | Default | Description |
| :---               | :---            | :---    | :---        |
| `logger`           | `SILENT`, `DISCRETE`, `WARNING`, `INFO`, `DEBUG`, `DEBUG2`, `DEBUG3` | `INFO` | Verbosity level of the logger |
| `progress_norm`    | `L1`, `L2`, `INF` | `L1` | Norm used for the progress measures |
| `residual_norm`    | `L1`, `L2`, `INF` | `INF` | Norm used for the residuals |

### Numerical options
| Option             | Type | Default | Description |
| :---               | :--- | :---    | :---        |
| `residual_scaling_threshold` | double | 100.0 | Scaling factor in stationarity and complementarity residuals |
| `protect_actual_reduction_against_roundoff` | bool | `false` | Whether the actual reduction is slightly modified to account for roundoff |
| `print_subproblem` | bool |  `false` | Whether the subproblem is printed in `DEBUG` mode |

## Globalization strategy options

| Option                     | Type   | Default | Description |
| :---                       | :---   | :---    | :---        |
| `armijo_decrease_fraction` | double | 1e-4    | Fraction of the predicted reduction that should be achieved by the actual reduction in the Armijo condition |
| `armijo_tolerance`         | double | 1e-9    | Minimum value of the predicted reduction in the Armijo condition, 0 otherwise |

## Switching method options

| Option                             | Type   | Default | Description |
| :---                               | :---   | :---    | :---        |
| `switching_delta`                  | double | 0.999   | Fraction of the constraint violation that should be achieved by the predicted reduction in the switching condition |
| `switching_infeasibility_exponent` | double | 2       | Exponent of the constraint violation in the switching condition |

## Filter method options

### String options

| Option        | Possible values           | Default    | Description |
| :---          | :---                      | :---       | :---        |
| `filter_type` | `standard`, `nonmonotone` | `standard` | Type of the filter data structure |

### Numerical options

| Option                                            | Type    | Default | Description |
| :---                                              | :---    | :---    | :---        |
| `filter_beta`                                     | double  | 0.999   | Fraction in the infeasibility sufficient reduction condition |
| `filter_gamma`                                    | double  | 0.001   | Slope in the objective sufficient reduction condition |
| `filter_ubd`                                      | double  | 1e2     | Minimum value for the initial upper bound on the infeasibility |
| `filter_fact`                                     | double  | 1.25    | Multiple of the initial infeasibility for the initial upper bound on the infeasibility |
| `filter_capacity`                                 | integer | 50      | Maximum number of filter entries |
| `filter_sufficient_infeasibility_decrease_factor` | double  | 0.9     | Infeasibility decrease factor in the infeasibility sufficient decrease condition |
| `nonmonotone_filter_number_dominated_entries`     | integer | 3       | Number of dominated filter entries (nonmonotone filter) |

## Funnel options

| Option                                          | Type    | Default | Description |
| :---                                            | :---    | :---    | :---        |
| `funnel_kappa`                                  | double  | 0.5     | Convex combination coefficient in funnel update rule |
| `funnel_beta`                                   | double  | 0.9999  | Fraction in the infeasibility sufficient reduction condition |
| `funnel_gamma`                                  | double  | 0.001   | Slope in the objective sufficient reduction condition |
| `funnel_ubd`                                    | double  | 1.0     | Minimum value for the initial upper bound on the infeasibility |
| `funnel_fact`                                   | double  | 1.5     | Multiple of the initial infeasibility for the initial upper bound on the infeasibility |
| `funnel_update_strategy`                        | integer | 1       | Rule for the funnel update (1, 2, or 3) |
| `funnel_require_acceptance_wrt_current_iterate` | bool    | `false` | Whether the trial iterate should improve upon the current iterate |

## Line search options

| Option                            | Type   | Default | Description |
| :---                              | :---   | :---    | :---        |
| `LS_backtracking_ratio`           | double | 0.5     | Decrease ratio of the step length for backtracking line search |
| `LS_min_step_length`              | double | 1e-12   | Minimum acceptable step length before failure is reported |
| `LS_scale_duals_with_step_length` | bool   | `true`  | Whether the Lagrange multipliers are scaled with the step length |

## Inertia correction options

| Option                                       | Type    | Default | Description |
| :---                                         | :---    | :---    | :---        |
| `regularization_failure_threshold`           | double  | 1e40    | Threshold for the primal inertia correction coefficient above which failure is reported |
| `primal_regularization_initial_factor`       | double  | 1e-4    | Initial value of primal inertia correction coefficient |
| `regularization_increase_factor`             | double  | 2       | Increase factor for the primal inertia correction coefficient |
| `dual_regularization_fraction`               | double  | 1e-8    | Fraction of the dual inertia correction parameter |
| `primal_regularization_lb`                   | double  | 1e-20   | Minimum value of the primal inertia correction coefficient upon decrease |
| `primal_regularization_decrease_factor`      | double  | 3.0     | Decrease factor for the primal inertia correction coefficient |
| `primal_regularization_fast_increase_factor` | double  | 100.0   | Fast increase factor for the primal inertia correction coefficient |
| `primal_regularization_slow_increase_factor` | double  | 8.0     | Slow increase factor for the primal inertia correction coefficient |
| `threshold_unsuccessful_attempts`            | integer | 8       | Number of unsuccessful attempts until inertia correction becomes more aggressive |

## Trust region options

| Option                          | Type   | Default | Description |
| :---                            | :---   | :---    | :---        |
| `TR_radius`                     | double | 10.0    | Initial value of the radius |
| `TR_increase_factor`            | double | 2       | Increase factor of the radius for successful iterations |
| `TR_decrease_factor`            | double | 2       | Decrease factor of the radius for unsuccessful iterations |
| `TR_aggressive_decrease_factor` | double | 4       | Decrease factor of the radius when errors occur |
| `TR_activity_tolerance`         | double | 1e-6    | Tolerance with which the trust-region constraint is considered active |
| `TR_min_radius`                 | double | 1e-7    | Minimum radius acceptable before failure is reported |
| `TR_radius_reset_threshold`     | double | 1e-4    | Smallest value to which the radius is reset at the beginning of an iteration |

## Feasibility restoration options

| Option                                                 | Type   | Default | Description |
| :---                                                   | :---   | :---    | :---        |
| `switch_to_optimality_requires_linearized_feasibility` | bool   | `true`  | Whether the switch to optimality phase requires the linearized constraints to be consistent |
| `l1_constraint_violation_coefficient`                  | double | 1       | Coefficient of the constraint violation in the $\ell_1$ relaxed problem |

## Barrier subproblem options

### String options

| Option             | Possible values | Default    | Description |
| :---               | :---            | :---       | :---        |
| `barrier_function` | `log`           | `log`      | Type of the barrier function |

### Numerical options

| Option                                 | Type   | Default | Description |
| :---                                   | :---   | :---    | :---        |
| `barrier_initial_parameter`            | double | 0.1     | Initial value of the barrier parameter |
| `barrier_default_multiplier`           | double | 1       | Initial value of the bound multipliers |
| `barrier_tau_min`                      | double | 0.99    | Coefficient of the fraction-to-boundary rule |
| `barrier_k_sigma`                      | double | 1e10    | Safeguard parameter for rescaling the bound multipliers |
| `barrier_k_mu`                         | double | 0.2     | Coefficient for multiplicative update of the barrier parameter |
| `barrier_theta_mu`                     | double | 1.5     | Coefficient for geometric update of the barrier parameter |
| `barrier_k_epsilon`                    | double | 10      | Scaling factor of the barrier parameter in the barrier update rule |
| `barrier_update_fraction`              | double | 10      | Fraction in the barrier update rule |
| `barrier_regularization_exponent`      | double | 0.25    | Exponent of the barrier parameter in the dual inertia correction parameter |
| `barrier_small_direction_factor`       | double | 10.0    | Multiple of the machine epsilon in the small step condition |
| `barrier_push_variable_to_interior_k1` | double | 1e-2    | Coefficient for the perturbation of the initial bounds |
| `barrier_push_variable_to_interior_k2` | double | 1e-2    | Coefficient for the perturbation of the initial bounds |
| `barrier_damping_factor`               | double | 1e-5    | Damping coefficient for single bounds |
| `least_square_multiplier_max_norm`     | double | 1e3     | Maximum accepted norm of the least-square multipliers |

## BQPD options

| Option                | Possible values         | Default     | Description |
| :---                  | :---                    | :---        | :---        |
| `BQPD_kmax_heuristic` | `filtersqp`, `minotaur` | `filtersqp` | Heuristic used to pick upper bound on nullspace size (`kmax`) |