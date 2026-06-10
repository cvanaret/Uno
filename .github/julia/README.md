# Julia tests

The Julia workflows run different solver configurations on MINLPTests.jl and MathOptInterface.jl instances.
We have partitioned the instances in order to target particular tests (e.g., inequality-constrained problems).

| Option                                   | unconstrained | bound-constrained | equality-constrained | general | LP | NLP |
|:-----------------------------------------|:--------------|:------------------|:---------------------|:--------|:---|:----|
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| ``                                       |               |                   |                      |         |    |     |
| `runtests_uno_pure_sqp_mumps`            |               |                   | x                    |         |    | x   |
| `runtests_uno_unconstrained_lbfgs_mumps` | x             |                   |                      |         |    | x   |
| `runtests_uno_unconstrained_mumps`       | x             |                   |                      |         | x  | x   |
