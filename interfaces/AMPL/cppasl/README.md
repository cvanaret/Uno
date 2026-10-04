# cppasl — a C++17 reader and evaluator for AMPL .nl files (ASL2 re-implementation)

cppasl reads AMPL `.nl` files (text `g` and binary `b`/`z`/`h`) and provides the ASL2 functionality needed by
nonlinear solvers — values, gradients, sparse Jacobian, sparse Lagrangian Hessian, Hessian-vector products, linear
operators, and `.sol` files — with the ASL2 thread-safety model (immutable model + one `EvaluationWorkspace` per
thread). On top of ASL2, it detects the algebraic structure of the problem exactly:

* every polynomial part of degree ≤ 2 is extracted at read time and stored as a lower-triangular COO Hessian, so that
  quadratic functions are evaluated with sparse matrix-vector products (no AD at all for LP/QP/QCQP);
* the problem is classified as LP, QP, QCQP or NLP (`MI` prefix with integers), and QP/QCQP convexity is decided
  (cheap certificates + LAPACK `dpotrf` per connected component).

## Usage

```cpp
#include "cppasl/nl_model.hpp"

const cppasl::NlModel model = cppasl::NlModel::read_file("problem.nl");
std::cout << model.classify().to_string() << '\n';              // e.g. "convex QCQP"
cppasl::EvaluationWorkspace workspace(model);                    // one per thread
const double f = model.evaluate_objective(workspace, x);
model.evaluate_objective_gradient(workspace, x, gradient);
model.evaluate_constraints(workspace, x, c);
model.evaluate_jacobian(workspace, x, jacobian_values);          // CSR: jacobian_row_starts(), jacobian_column_indices()
model.evaluate_lagrangian_hessian(workspace, x, sigma, y, h);    // lower CSC: hessian_column_starts(), hessian_row_indices()
model.evaluate_lagrangian_hessian_vector_product(workspace, x, sigma, y, v, hv);   // matrix-free
model.multiply_jacobian(jacobian_values, v, Jv);                 // J v, J^T v, W v on evaluated values
model.objective_quadratic_form();                                // H of the quadratic part (lower COO)
model.write_solution("problem.sol", "message", x, y, 0);
```

Conventions: 0-based indices, Lagrangian `sigma ∇²f + Σ y_i ∇²c_i` (objective selected by
`ReaderOptions::hessian_objective_index`), infinite bounds are `±infinity`.

## Design

* **Streaming reader** (`src/nl_reader.cpp`): the file is memory-mapped; one templated segment reader per encoding
  (text: hand-written integer parser and `std::from_chars` for reals; binary: `memcpy`). Expressions are parsed with
  an explicit stack (arbitrarily deep expressions cannot overflow the call stack) into a flat arena.
* **Decomposition at read time** (`src/function_decomposition.cpp`, `src/model_builder.cpp`): each function is split,
  through top-level sums/negations/scalings and defined variables, into constant + linear + ½xᵀHx + Σ scale·element.
  Degree analysis is memoized; products of affine forms and squares are expanded symbolically; duplicates are merged
  with a two-pass counting sort. Right after a function is processed, its expression is discarded (only defined
  variables stay resident), so memory is proportional to the output, not to the file.
* **AD tapes** (`src/tape.cpp`): the non-quadratic elements are compiled into one flat, topologically ordered
  instruction array with global operand indices. Constant operands are folded into unary instructions (`u+c`, `c*u`,
  `u^2`, `u^n` by multiplications, `sqrt`), `x ± c` is a single shifted variable leaf, operation shapes are table
  lookups.
  * values + first partials in one forward sweep (as ASL), for a whole function or a whole run of consecutive
    constraints at once;
  * gradients and Jacobians by **reverse programs**: per function, a list of edges `adj[t] += adj[s] * partial[p]`
    by decreasing source (ASL's `derp` lists, 12 bytes per edge instead of 24), branch-free; the first edge into an
    instruction assigns, so adjoints are never zeroed; dead instructions (under piecewise-constant operations) get no
    edges; variable leaves then scatter to their gradient/Jacobian positions;
  * element Hessians by vector forward-over-reverse (8 directions per sweep pair), packed column-major like the CSC;
  * group elements `phi(sum_k t_k)` (log-sum-exp, norms, powers of sums, ...) with operands compiled as
    self-contained instruction ranges: `H = phi'' g g^T + phi' sum_k H_k`, each `H_k` swept on its own range (ASL's
    group partial separability);
  * **shared defined variables** (AMPL common expressions used by several functions, i.e. not `c1`/`o1`): evaluated
    once per point by their own tapes (one sweep for all of them), with their gradients (one reverse program); the
    functions read them through `DefinedValue` leaves and use the chain rule: Jacobian `+= adj_v grad v`; Hessian
    `J^T H_phi J` per element plus `(sum of the multiplier-weighted adj_v) H_v` added once per defined variable;
    Hessian-vector products seed the defined variables from all functions, then one seeded pass over their tapes.
    Elements sharing single-function defined variables are merged (bounded), so these are evaluated once per function.
* **Precomputed scatter maps**: linear coefficients are aligned with the CSR Jacobian; every quadratic entry and
  element-variable pair knows its Jacobian position and its Hessian slot, so the Jacobian/Hessian evaluations are
  pure scatter loops (e.g. the Hessian of a QCQP is `h[slot[q]] += w * H[q]`).
* **Classification** (`src/problem_classification.cpp`): `min f` needs PSD H (`max`: NSD); `g ≤ u` PSD, `g ≥ l` NSD,
  two-sided or equality with H ≠ 0 nonconvex. Indefiniteness certificates: negative diagonal, 2×2 minor; PSD
  certificate: diagonal dominance; otherwise dense Cholesky (`dpotrf`) of `H + δI`, `δ = 1e-8·max(1, max|H_ii|)`, per
  connected component of the sparsity graph. Components larger than `maximum_dense_dimension` (4000) give `Unknown`.

## Tests

`build/cppasl_tests [filter]` (16 test cases, ~800 checks): hand-written text file (HS071), text/binary equality,
bounds/suffixes/complementarity/integrality, `.sol` output, 11 AMPL encodings of quadratics compared with the AD path
(`detect_quadratic_structure = false`), all unary/binary/n-ary operations with finite-difference checks of gradients,
Jacobians, Hessians and HVPs, defined variables, workspace caching, classification (hidden indefiniteness,
maximization, one/two-sided constraints, 40000-variable block structure), and large instances (QCQP up to n = 100000,
86 MB, read in 0.6 s; partially separable NLP with n = 100000). `CPPASL_TEST_OUTPUT=<dir>` saves the test `.nl` files
for `compare_with_asl2`.

## Comparison with ASL2 (`compare_with_asl2`, `benchmarks/run_benchmarks.sh`)

Same machine (1 core), ASL2 = `pfgh_read` + `sphsetup` (upper triangle), best of 6 runs, alternating points.
Instances from Mittelmann's `qcqp.mod` (generator `generate_qcqp_nl`; AMPL's RNG is not reproduced).
All values agree entry by entry (max relative difference ≤ 1.4e-12), with identical Jacobian and Hessian patterns.

| instance (n, m, file) | read (+ ASL `sphsetup`) | f | ∇f | c | J | ∇²L | ∇²L·v | peak RSS |
|---|---|---|---|---|---|---|---|---|
| A: 500, 120, 82 MB, convex QCQP | 2.80 s → **0.39 s** (7.1×, 8.7× with setup) | 26× | 18× | 70× | 15× | 520× | 50× | 942 → 182 MB |
| B: 1500, 10508, 48 MB, convex QCQP | 1.48 s → **0.25 s** (6.0×, 7.7×) | 27× | 10× | 43× | 15× | 260× | 100× | 581 → 116 MB |
| 750, 120, 190 MB, sd = 0, nonconvex QCQP | 6.97 s → **1.05 s** (6.6×, 8.0×) | 26× | 18× | 67× | 15× | 564× | 40× | 2172 → 424 MB |

(cppasl's RSS includes the memory-mapped file.) Classification takes 0.14–0.34 s (dense `dpotrf` on the 1500×1500
blocks), 1 ms when a 2×2 certificate proves nonconvexity.

### Large nonlinear instances (`generate_nlp_nl`)

Speedup of cppasl over ASL2 per call (≥ 1 means faster), values agreeing to ≤ 6e-14 (except ASL's own errors below):

| instance (n = 500,000 unless stated) | read¹ | f | ∇f | c | J | ∇²L | ∇²L·v |
|---|---|---|---|---|---|---|---|
| chain: 1.5M tiny elements, 500k constraints | 1.4× | 1.2× | 1.4× | 1.2× | 1.0× | 1.4× | 1.5× |
| defined: defined variables in the objective and one constraint | 1.4× | 1.4× | 1.9× | 1.0× | 1.1× | 2.5× | 1.9× |
| shared (n = 100,000): each defined variable (20 vars) in 20 constraints | 2.5× | 1.0× | 1.6× | 2.0× | 4.0× | 10× | 4.2ײ |
| blocks: log-sum-exp + norm groups of 20 variables | 2.1× | 1.8× | 1.8× | 2.3× | 2.1× | 2.2× | 2.2× |
| dense: one element over all n = 3000 variables | 1.1× | 2.2× | 2.0× | 1.4× | 1.9× | 1.2× | 1.1× |

¹ against ASL's read + `sphsetup` (cppasl builds the Hessian structure while reading).
² ASL2's `hvcomp` is wrong on this model (see below); timing comparison only.

### Memory layout

All graphs are flat arrays (no linked structures, unlike ASL's `cgrad`/`ograd` lists and pointer-linked expressions):
the expression arena (nodes + 32-bit operand indices), one contiguous instruction array for all tapes (elements of a
function and runs of constraints adjacent, swept sequentially), values/partials/adjoints indexed like the
instructions, K-direction tangent blocks interleaved per instruction, CSR/CSC matrices and 32-bit scatter slots
(bucketed by counting sorts). The remaining indirect accesses are the scatter into the Hessian slots and the
operand lists of n-ary sums.

### Discrepancies found in ASL2 (cppasl verified by finite differences)

* `asinh` (o50), second-order code (`eval2.c`, `OPasinh_g` and `OP_asinh1`): the value `log(t + sqrt(t^2+1))` is
  computed before taking `|t|` and then negated for negative arguments, so the value used by Jacobians and Hessians
  has the wrong sign for t < 0.
* `hvcomp` with defined variables shared by many constraints (the `shared` instance): ASL2's Hessian-vector product
  differs from ASL2's own `sphes` Hessian times the vector by up to 38% (cppasl's HVP matches its Hessian to 2e-16,
  and its Hessian matches `sphes`).

* `signpow` (o80): ASL2's Hessian code aborts (`Bad *o = 202 in hv_back`).
* `x mod y` (o4): ASL2's second derivatives are zero; the mixed derivative of `rem(10 x0, 3) * x1` is 10.
* piecewise-linear terms (o64): ASL assumes a variable argument (the Hessian of `<<...>> (x0*x1)` is missing), and for
  a negative argument inside the segment containing 0 (e.g. breakpoints -0.5, 0.5 and x = -0.3) it uses the slope of
  the segment to the left (ASL1 has the same code; this is harmless if AMPL always emits a breakpoint at 0).

## Limitations

Imported functions (`F`), logical constraints (`L`), string arguments and counting/logical operators outside of
conditions are rejected with an error; binary files must have the native byte order; element Hessian patterns are
dense over the element's variables (ASL can be sparser, e.g. `x0/x1` or `min`); convexity of connected quadratic forms
with more than 4000 variables that are not diagonally dominant is reported as `Unknown` (a sparse Cholesky would lift
this); the shared defined variables are all evaluated as soon as one function needs one of them (ASL tracks which
common expressions each function needs).

## Build

```
cmake -B build -DCMAKE_BUILD_TYPE=Release [-DCPPASL_USE_LAPACK=ON] [-DCPPASL_ASL2_ROOT=/path/to/asl]
cmake --build build && build/cppasl_tests
benchmarks/run_benchmarks.sh build /tmp/instances    # needs CPPASL_ASL2_ROOT with build/lib/libasl2.a
```
