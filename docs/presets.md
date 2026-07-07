# Presets

Uno implements **presets**, that is combinations of ingredients and hyperparameters that mimic existing solvers. Uno currently supports:

* `filtersqp` mimics filterSQP (trust-region feasibility restoration filter SQP method with exact Hessian);
* `ipopt` mimics IPOPT (line-search feasibility restoration filter barrier method with exact Hessian and primal-dual inertia correction).

The values listed below are taken from `uno/options/Presets.cpp`.

The **default preset** `auto` decides between `filtersqp` and `ipopt`, depending on the properties of the problem. The following (simple) oracle is currently used:

- if $2000 \le n + m$ or $50000 \le nnz(\nabla c) + nnz(\nabla^2 \mathcal{L})$, pick the `ipopt` preset
- otherwise, pick the `filtersqp` preset

## `filtersqp` preset

| Option                                                 | Type   | Value                     |
|:-------------------------------------------------------|:-------|:--------------------------|
| `constraint_relaxation_strategy`                       | string | `feasibility_restoration` |
| `inequality_handling_method`                           | string | `inequality_constrained`  |
| `hessian_model`                                        | string | `exact`                   |
| `inertia_correction_strategy`                          | string | `none`                    |
| `globalization_mechanism`                              | string | `TR`                      |
| `globalization_strategy`                               | string | `fletcher_filter_method`  |
| `filter_type`                                          | string | `standard`                |
| `progress_norm`                                        | string | `L1`                      |
| `residual_norm`                                        | string | `L2`                      |
| `TR_radius`                                            | double | 10                        |
| `l1_constraint_violation_coefficient`                  | double | 1                         |
| `primal_tolerance`                                     | double | 1e-6                      |
| `dual_tolerance`                                       | double | 1e-6                      |
| `switch_to_optimality_requires_linearized_feasibility` | bool   | true                      |
| `protect_actual_reduction_against_roundoff`            | bool   | false                     |

If not provided, the QP solver is chosen automatically from the available QP solvers (if any).

## `ipopt` preset

| Option                                                 | Type   | Value                     |
|:-------------------------------------------------------|:-------|:--------------------------|
| `constraint_relaxation_strategy`                       | string | `feasibility_restoration` |
| `inequality_handling_method`                           | string | `interior_point`          |
| `barrier_function`                                     | string | `log`                     |
| `hessian_model`                                        | string | `exact`                   |
| `inertia_correction_strategy`                          | string | `primal_dual`             |
| `globalization_mechanism`                              | string | `LS`                      |
| `globalization_strategy`                               | string | `waechter_filter_method`  |
| `filter_type`                                          | string | `standard`                |
| `filter_beta`                                          | double | 0.99999                   |
| `filter_gamma`                                         | double | 1e-8                      |
| `switching_delta`                                      | double | 1                         |
| `filter_ubd`                                           | double | 1e4                       |
| `filter_fact`                                          | double | 1e4                       |
| `filter_switching_infeasibility_exponent`              | double | 1.1                       |
| `armijo_decrease_fraction`                             | double | 1e-8                      |
| `LS_backtracking_ratio`                                | double | 0.5                       |
| `LS_min_step_length`                                   | double | 5e-7                      |
| `barrier_tau_min`                                      | double | 0.99                      |
| `barrier_damping_factor`                               | double | 1e-5                      |
| `l1_constraint_violation_coefficient`                  | double | 1000                      |
| `progress_norm`                                        | string | `L1`                      |
| `residual_norm`                                        | string | `INF`                     |
| `primal_tolerance`                                     | double | 1e-8                      |
| `dual_tolerance`                                       | double | 1e-8                      |
| `loose_primal_tolerance`                               | double | 1e-6                      |
| `loose_dual_tolerance`                                 | double | 1e-6                      |
| `loose_tolerance_iteration_threshold`                  | int    | 15                        |
| `switch_to_optimality_requires_linearized_feasibility` | bool   | false                     |
| `LS_scale_duals_with_step_length`                      | bool   | true                      |
| `protect_actual_reduction_against_roundoff`            | bool   | true                      |

If not provided, the linear solver is chosen automatically from the available linear solvers (if any).