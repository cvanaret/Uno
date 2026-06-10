# Julia tests

The Julia workflows run different solver configurations on MINLPTests.jl and MathOptInterface.jl instances.
We have partitioned the instances in order to target particular tests (e.g., inequality-constrained problems).

| Script                                          | unconstrained | bound-constrained | equality-constrained | general | LP | NLP |
|:------------------------------------------------|:--------------|:------------------|:---------------------|:--------|:---|:----|
| `runtests_uno_filterslp_bqpd`                   |               |                   |                      |         |    |     |
| `runtests_uno_filterslp_highs`                  |               |                   |                      |         |    |     |
| `runtests_uno_filtersqp_bqpd`                   |               |                   |                      |         |    |     |
| `runtests_uno_filtersqp_bqpd_identity`          |               |                   |                      |         |    |     |
| `runtests_uno_filtersqp_bqpd_ls_lbfgs`          |               |                   |                      |         |    |     |
| `runtests_uno_filtersqp_bqpd_lsr1`              |               |                   |                      |         |    |     |
| `runtests_uno_filtersqp_bqpd_regularization_ls` |               |                   |                      |         |    |     |
| `runtests_uno_filtersqp_bqpd_regularization_tr` |               |                   |                      |         |    |     |
| `runtests_uno_ipopt_ma27`                       |               | x                 |                      | x       | x  | x   |
| `runtests_uno_ipopt_ma57`                       |               | x                 |                      | x       | x  | x   |
| `runtests_uno_ipopt_mumps`                      |               | x                 |                      | x       | x  | x   |
| `runtests_uno_ipopt_ssids`                      |               | x                 |                      | x       | x  | x   |
| `runtests_uno_ipopt_mumps_identity`             |               | x                 |                      | x       |    | x   |
| `runtests_uno_ipopt_mumps_lbfgs`                |               | x                 |                      | x       |    | x   |
| `runtests_uno_pure_sqp_mumps`                   |               |                   | x                    |         |    | x   |
| `runtests_uno_unconstrained_lbfgs_mumps`        | x             |                   |                      |         |    | x   |
| `runtests_uno_unconstrained_mumps`              | x             |                   |                      |         | x  | x   |

`runtests_uno_filtersqp_highs` runs only on convex `MINLPTests.nlp_cvx_expr`.