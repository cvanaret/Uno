# SciPy interface

`unopy.scipy.minimize` solves an optimization problem with Uno through the API of
[`scipy.optimize.minimize`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html):
the objective, its derivatives, the bounds and the constraints are given as in SciPy,
and the result is a SciPy `OptimizeResult`.
Code written for `scipy.optimize.minimize` therefore runs with Uno by changing the import and the method.

It is part of [unopy](python.md) and its dependencies are installed with the `scipy` extra:
```bash
pip install "unopy[scipy]"
```

## Solving a problem

The methods are the [Uno presets](../presets.md): `"filtersqp"` (the default), `"filterslp"`, `"funnelsqp"` and `"ipopt"`.

```python
import numpy as np
from scipy.optimize import rosen, rosen_der, rosen_hess
from unopy.scipy import minimize

x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
res = minimize(rosen, x0, jac=rosen_der, hess=rosen_hess, method="filtersqp")
print(res.x)        # [1. 1. 1. 1. 1.]
print(res.success)  # True
print(res.message)  # Optimization status: SUCCESS; solution status: FEASIBLE_KKT_POINT.
```

The arguments are those of `scipy.optimize.minimize`:
`fun`, `x0`, `args`, `method`, `jac`, `hess`, `bounds`, `constraints`, `tol`, `callback` and `options`.

## Derivatives

The gradient `jac` can be:

- a callable `jac(x, *args)` returning the gradient;
- `True`, if `fun` returns both the objective and its gradient;
- a finite-difference scheme: `"2-point"` (the default when `jac` is not given), `"3-point"` or `"cs"` (complex step).

The Hessian `hess` can be:

- a callable `hess(x, *args)` returning the Hessian of the objective as a dense array, a `scipy.sparse` matrix or a `LinearOperator`.
  The Hessian of each `NonlinearConstraint` is then read from its `hess` attribute,
  `hess(x, v)`, returning the sum of the Hessians of its components weighted by the multipliers `v`.
  Uno only uses Hessian-vector products, so the Lagrangian Hessian is never built as a dense matrix;
- a `scipy.optimize.BFGS()` or `SR1()` instance: Uno then uses its limited-memory quasi-Newton model,
  L-BFGS or L-SR1 (the option `quasi_newton_memory_size` sets the memory);
- `None` (the default): Uno uses its L-BFGS model.

A callable `hess` is used by `"filtersqp"` and `"funnelsqp"`.
It is ignored by `"filterslp"`, which uses no Hessian, and when the option `hessian_model` is not `"exact"`.
It is not supported by `"ipopt"`, which needs an explicit Hessian matrix: use `BFGS()` with this method.

## Bounds and constraints

The bounds are a `scipy.optimize.Bounds` instance or a sequence of `(min, max)` pairs, where `None` means no bound.
The constraints are `NonlinearConstraint` and `LinearConstraint` instances or dictionaries,
alone or in a list.

```python
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint
from scipy.optimize import rosen, rosen_der, rosen_hess
from unopy.scipy import minimize

bounds = Bounds([0.0, -0.5], [1.0, 2.0])
linear_constraint = LinearConstraint([[1.0, 2.0], [2.0, 1.0]], [-np.inf, 1.0], [1.0, 1.0])
nonlinear_constraint = NonlinearConstraint(
    lambda x: [x[0] ** 2 + x[1], x[0] ** 2 - x[1]],
    -np.inf,
    1.0,
    jac=lambda x: [[2 * x[0], 1.0], [2 * x[0], -1.0]],
    hess=lambda x, v: 2 * (v[0] + v[1]) * np.array([[1.0, 0.0], [0.0, 0.0]]),
)
res = minimize(
    rosen,
    np.array([0.5, 0.0]),
    jac=rosen_der,
    hess=rosen_hess,
    method="funnelsqp",
    bounds=bounds,
    constraints=[linear_constraint, nonlinear_constraint],
)
print(res.x)  # [0.41494432 0.17011137]
```

A similar problem with dictionary constraints, whose Jacobians are approximated by finite differences,
solved with the IPOPT preset:

```python
constraints = [
    {"type": "ineq", "fun": lambda x: np.array([1 - x[0] - 2 * x[1], 1 - x[0] ** 2 - x[1]])},
    {"type": "eq", "fun": lambda x: np.array([2 * x[0] + x[1] - 1])},
]
res = minimize(
    rosen,
    np.array([0.5, 0.0]),
    method="ipopt",
    bounds=[(0, 1), (-0.5, 2.0)],
    constraints=constraints,
    tol=1e-8,
)
print(res.x)  # [0.41494432 0.17011137]
```

A dictionary constraint `{"type": "ineq", "fun": g}` means `g(x) >= 0`, as in SciPy;
its `"jac"` entry is optional, and its `"args"` entry gives extra arguments to `fun` and `jac`.

### Sparse Jacobians

The Jacobian of a constraint, and the matrix `A` of a `LinearConstraint`,
can be a `scipy.sparse` matrix: it is passed to Uno without being converted to a dense array.
The sparsity pattern is fixed for the whole optimization.
It is the pattern of the Jacobian at `x0`, which is evaluated before the optimization,
or the nonzero entries of a dense `A`.
A sparse Jacobian must therefore store every entry that can be nonzero during the optimization,
even if it is zero at `x0`; for instance, a sparse matrix built from a dense array drops its zeros.
A COO matrix built from its entries keeps them:

```python
import numpy as np
import scipy.sparse as sp
from scipy.optimize import NonlinearConstraint, rosen, rosen_der
from unopy.scipy import minimize

n = 100
rows = np.arange(n - 1)

def constraint_jacobian(x):
    # The entries 2 * x[i] are stored even when they are zero, e.g. at x0 = 0.
    return sp.coo_array(
        (
            np.concatenate([2 * x[:-1], np.ones(n - 1)]),
            (np.concatenate([rows, rows]), np.concatenate([rows, rows + 1])),
        ),
        shape=(n - 1, n),
    )

constraint = NonlinearConstraint(
    lambda x: x[:-1] ** 2 + x[1:], -np.inf, 2.0, jac=constraint_jacobian
)
res = minimize(rosen, np.zeros(n), jac=rosen_der, method="ipopt", constraints=constraint)
```

A nonzero entry outside the pattern raises a `ValueError`.

## Options

`tol` sets both the `primal_tolerance` and `dual_tolerance` options of Uno.
The other [Uno options](../options.md) are passed in the `options` dictionary,
with `"maxiter"` and `"max_iter"` as aliases of `"max_iterations"`:

```python
res = minimize(
    rosen,
    x0,
    jac=rosen_der,
    tol=1e-10,
    options={"maxiter": 500, "logger": "SILENT", "TR_radius": 0.5},
)
```

An unknown option or a value of the wrong type raises an error before the optimization.

## Callback

The callback is called after each iteration.
With a parameter named `intermediate_result`, it is passed an `OptimizeResult` with the current point `x` and objective `fun`;
otherwise, it is passed a copy of the current point.
Raising `StopIteration` in the callback stops Uno cleanly:

```python
def callback(intermediate_result):
    print(intermediate_result.fun)
    if intermediate_result.fun < 1e-2:
        raise StopIteration

res = minimize(rosen, x0, jac=rosen_der, callback=callback)
print(res.success, res.status)  # False 99
print(res.message)              # `callback` raised `StopIteration`.
```

## Errors in the user functions

An exception raised by `fun`, `jac`, `hess`, a constraint function or the callback (other than `StopIteration`)
stops Uno, and the same exception is raised by `minimize`:

```python
def fun(x):
    if x[0] > 1.2:
        raise RuntimeError(f"cannot evaluate at {x}")
    return rosen(x)

minimize(fun, np.full(5, 0.5), jac=rosen_der)  # RuntimeError: cannot evaluate at [...]
```

The user functions are no longer called once one of them has raised an exception.

## Result

`minimize` returns an `OptimizeResult` with the attributes of SciPy:

| Attribute | Description |
|-----------|-------------|
| `x` | the solution |
| `fun` | the objective value at `x` |
| `success` | whether Uno returned the optimization status `SUCCESS` with a feasible solution (`FEASIBLE_KKT_POINT`, `FEASIBLE_FJ_POINT` or `FEASIBLE_SMALL_STEP`) |
| `status` | the optimization status as an int: 0 when Uno terminated normally (`success` can still be False at an infeasible point), 99 when the callback raised `StopIteration` |
| `message` | the optimization status and the solution status |
| `nit` | the number of iterations |
| `nfev` | the number of calls to `fun`, finite differences included |
| `njev` | the number of evaluations of the gradient |
| `nhev` | the number of calls to `hess` |
| `maxcv` | the maximum violation of the bounds and the constraints at `x` |

and the following ones:

| Attribute | Description |
|-----------|-------------|
| `optimization_status`, `solution_status` | the Uno statuses, as `unopy.OptimizationStatus` and `unopy.SolutionStatus` |
| `constr` | the values of the constraints at `x`, one array per constraint |
| `v` | the Lagrange multipliers of the constraints, one array per constraint, followed by those of the bounds if `bounds` is given |
| `lower_bound_multipliers`, `upper_bound_multipliers` | the multipliers of the lower and upper bounds of the variables |

The multipliers follow the convention of SciPy's `trust-constr` method:
the Lagrangian is $L(x, v, z) = f(x) + \sum_i v_i^T c_i(x) + z^T x$ and its gradient is zero at a KKT point.
A positive multiplier means that the upper bound of the constraint is active and a negative one that its lower bound is active.
