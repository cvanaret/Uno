## Iterative refinement

Iterative refinement is performed after each solve with a direct linear solver. The residuals are computed on the uncondensed (possibly unsymmetric) Newton system, i.e. before the condensation of its condensed blocks (e.g. the bound duals of the primal-dual interior-point method); without condensation, the uncondensed system is the augmented system. The residuals use the factorized matrix (regularization included) without the Schur complement of the condensed blocks. Iterative refinement is not performed with limited-memory quasi-Newton Hessians (L-BFGS, L-SR1), since the augmented matrix then does not contain the full Hessian.

If refinement fails, the linear system scaling is switched on (if available and not yet on, see [Linear system scaling](#linear-system-scaling)); otherwise, the pivoting tolerance of the linear solver (MA27, MA57, MUMPS) is increased (at most once per matrix) and the matrix is refactorized; its inertia is checked again and the regularization increased if needed. The increase of the pivoting tolerance is permanent: all subsequent factorizations use the stricter tolerance, up to a solver-specific maximum (1e-4 for MA27 and MA57, 0.1 for MUMPS, as in IPOPT). If the tolerance cannot be increased and the residual ratio is above `iterative_refinement_residual_ratio_singular`, the matrix is deemed singular and the inertia correction strategy increases the regularization. Second-order corrections never modify the factorization (their refined solution is accepted as is).

Setting `iterative_refinement_max_steps` to 0 disables the whole mechanism (no residual computation, no pivoting tolerance increase, no singularity heuristic). Unlike IPOPT's `max_refinement_steps` (which allows up to `max_refinement_steps` + 1 steps), `iterative_refinement_max_steps` is the exact maximum number of refinement steps.

## Linear system scaling

The matrices of the direct linear solver can be scaled symmetrically (S K S, S a positive diagonal matrix: the inertia is preserved), independently of the linear solver, with `linear_system_scaling` (`none` by default). It replaces the solvers' own scalings: `MA57_use_scaling` and `MUMPS_use_scaling` are then ignored. Computing a scaling at every factorization is expensive, so by default the scaling is switched on only when the factorizations seem inaccurate (as IPOPT's `linear_scaling_on_demand`), then kept until the end:
- iterative refinement failed (before the pivoting tolerance is increased);
- optionally, the first factorization of a matrix delayed more than `linear_system_scaling_delayed_pivots_threshold` × dimension pivots (MA57, MA86, MUMPS, SSIDS); the scaling then starts at the next matrix. This trigger is disabled by default: the delayed pivots of KKT matrices are mostly structural (2×2 pivots around the zero (2, 2) block) and do not indicate inaccurate factorizations;
- the inertia correction failed (last resort, before the inertia correction gives up): the inertia correction restarts once with scaling.

The scaling factors are computed at most once per matrix: the refactorizations of the same matrix (inertia correction retries, pivoting tolerance increase) and the solves with an existing factorization (second-order corrections, iterative refinement) reuse them. Well-behaved problems never trigger the scaling and pay nothing.

`ruiz` is the symmetric Ruiz infinity-norm equilibration (Ruiz 2001; Knight, Ruiz & Uçar 2014).
