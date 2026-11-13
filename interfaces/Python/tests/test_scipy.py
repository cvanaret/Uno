# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""Tests of the scipy.optimize.minimize-like interface unopy.scipy."""

from __future__ import annotations

import inspect
import re
import warnings
from pathlib import Path

import pytest

pytest.importorskip("scipy")

import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import unopy  # noqa: E402
from scipy.optimize import BFGS  # noqa: E402
from scipy.optimize import Bounds  # noqa: E402
from scipy.optimize import HessianUpdateStrategy  # noqa: E402
from scipy.optimize import LinearConstraint  # noqa: E402
from scipy.optimize import NonlinearConstraint  # noqa: E402
from scipy.optimize import OptimizeResult  # noqa: E402
from scipy.optimize import SR1  # noqa: E402
from scipy.optimize import minimize as scipy_minimize  # noqa: E402
from scipy.optimize import rosen  # noqa: E402
from scipy.optimize import rosen_der  # noqa: E402
from scipy.optimize import rosen_hess  # noqa: E402
from scipy.sparse.linalg import LinearOperator  # noqa: E402

from unopy.scipy import minimize  # noqa: E402
from unopy.scipy._constraints import _JacobianBlock  # noqa: E402
from unopy.scipy._constraints import _normalize_constraints  # noqa: E402
from unopy.scipy._finite_differences import _approximate_jacobian  # noqa: E402
from unopy.scipy._guard import _EvaluationErrorGuard  # noqa: E402
from unopy.scipy._hessian import _LagrangianHessianOperator  # noqa: E402
from unopy.scipy._objective import _Objective  # noqa: E402


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp", "ipopt"])
def test_rosen(method):
    res = minimize(
        rosen,
        jac=rosen_der,
        x0=np.zeros(2),
        method=method,
        options={"max_iterations": 10000},
        tol=1e-4,
    )
    assert res.fun < 1e-4
    assert res.success


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp", "ipopt"])
def test_rosen_bnds(method):
    res = minimize(
        rosen,
        jac=rosen_der,
        x0=np.zeros(2),
        method=method,
        bounds=[[-1.0, 0.5], [-1, 1]],
        options={"max_iterations": 10000},
    )
    assert res.success


@pytest.mark.parametrize(
    "c_type", ["lambda", "NonLinearConstraint", "LinearConstraint"]
)
def test_rosen_constr(c_type):
    if c_type == "NonLinearConstraint":
        constr = NonlinearConstraint(
            fun=lambda x: 1.0 - 2 * x,
            jac=lambda x: -2 * np.eye(2),
            lb=np.zeros(2),
            ub=np.full(2, float("Inf")),
        )
    elif c_type == "LinearConstraint":
        constr = LinearConstraint(A=np.eye(2), ub=np.full(2, 0.5), lb=-float("Inf"))
    else:
        constr = {
            "type": "ineq",
            "fun": lambda x: 1.0 - 2 * x,
            "jac": lambda x: -2 * np.eye(2),
        }
    res = minimize(rosen, jac=rosen_der, x0=np.zeros(2), constraints=(constr,))

    assert res.success
    # Every constraint is x <= 0.5, active for x[0] at the solution.
    assert np.all(res.x <= 0.5 + 1e-6)
    assert np.allclose(res.x, [0.5, 0.25], atol=1e-4)

@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp", "ipopt"])
def test_rosen_constr2(method):
    res = minimize(
        rosen,
        np.array([1.3, 0.7, 0.8, 0.3]),
        jac=rosen_der,
        method=method,
        tol=1e-3,
        constraints=[
            NonlinearConstraint(
                fun=lambda x: x**2,
                jac=lambda x: 2 * np.diag(x),
                lb=np.full(4, 0.1),
                ub=np.full(4, 0.8),
            ),
        ],
        options={"max_iterations": 10000},
    )
    assert res.success

class UserError(Exception):
    """An error raised by a user function."""


def raise_after(function, number_of_calls):
    """Make a function raise a UserError after a number of calls.

    Args:
        function: The function to wrap.
        number_of_calls: The number of calls before raising.

    Returns:
        The wrapped function and the list of the points where it was called.
    """
    calls = []

    def wrapped_function(x, *args):
        calls.append(np.array(x))
        if len(calls) > number_of_calls:
            raise UserError(f"Failure at call {len(calls)}.")
        return function(x, *args)

    return wrapped_function, calls


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp", "ipopt"])
@pytest.mark.parametrize("failing", ["fun", "jac", "constraint", "constraint_jac"])
def test_user_exception(method, failing):
    """Check that an exception raised by a user function stops Uno and is re-raised."""
    functions = {
        "fun": rosen,
        "jac": rosen_der,
        "constraint": lambda x: x**2,
        "constraint_jac": lambda x: 2 * np.diag(x),
    }
    functions[failing], calls = raise_after(functions[failing], 3)
    with pytest.raises(UserError, match="Failure at call 4."):
        minimize(
            functions["fun"],
            np.array([1.3, 0.7, 0.8]),
            jac=functions["jac"],
            method=method,
            constraints=[
                NonlinearConstraint(
                    fun=functions["constraint"],
                    jac=functions["constraint_jac"],
                    lb=np.full(3, 0.1),
                    ub=np.full(3, 0.8),
                ),
            ],
            options={"max_iterations": 10000},
        )
    # The user function is no longer called once it has raised.
    assert len(calls) == 4


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp"])
def test_user_exception_hess(method):
    """Check that an exception raised by the Hessian stops Uno and is re-raised."""
    hess, calls = raise_after(rosen_hess, 2)
    with pytest.raises(UserError, match="Failure at call 3."):
        minimize(
            rosen,
            np.array([1.3, 0.7, 0.8]),
            jac=rosen_der,
            hess=hess,
            method=method,
            options={"max_iterations": 10000},
        )
    assert len(calls) == 3


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp"])
def test_rosen_constr2_hess(method):
    res = minimize(
        rosen,
        np.array([1.3, 0.7, 0.8]),
        jac=rosen_der,
        hess=rosen_hess,
        method=method,
        tol=1e-3,
        constraints=[
            NonlinearConstraint(
                fun=lambda x: x**2,
                jac=lambda x: 2 * np.diag(x),
                hess=lambda x, v: 2 * np.diag(v),
                lb=np.full(3, 0.1),
                ub=np.full(3, 0.8),
            ),
        ],
        options={"max_iterations": 10000},
    )
    assert res.success


def _sphere(x):
    return float(((x - 2.0) ** 2).sum())


def _sphere_der(x):
    return 2 * (x - 2.0)


def _sphere_hess(x):
    return 2 * np.eye(len(x))


# Two blocks with non-symmetric Jacobians, so that the order of the Jacobian values matters.
_A1 = np.array([[1.0, 2.0, 3.0]])
_A2 = np.array([[1.0, 0.0, -4.0], [0.5, 1.0, 1.0]])


def _blocks_constraints(with_hess=False):
    kwargs1 = {"hess": lambda x, v: 2 * v[0] * np.diag([0.0, 0.0, 1.0])} if with_hess else {}
    kwargs2 = {"hess": lambda x, v: np.zeros((3, 3))} if with_hess else {}
    return [
        NonlinearConstraint(
            fun=lambda x: _A1 @ x + x[2] ** 2,
            jac=lambda x: _A1 + np.array([[0.0, 0.0, 2 * x[2]]]),
            lb=-np.inf,
            ub=2.0,
            **kwargs1,
        ),
        NonlinearConstraint(
            fun=lambda x: _A2 @ x, jac=lambda x: _A2, lb=-1.0, ub=1.0, **kwargs2
        ),
    ]


def _reference_blocks_solution():
    return scipy_minimize(
        _sphere,
        np.zeros(3),
        jac=_sphere_der,
        method="SLSQP",
        constraints=_blocks_constraints(),
        options={"ftol": 1e-12},
    ).x


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp", "ipopt"])
def test_several_constraint_blocks(method):
    """Check a problem with several constraint blocks of different sizes.

    The bounds of the NonlinearConstraint are scalars broadcast to the constraint size.
    """
    res = minimize(
        _sphere,
        np.zeros(3),
        jac=_sphere_der,
        method=method,
        tol=1e-8,
        constraints=_blocks_constraints(),
        options={"max_iterations": 1000},
    )
    assert res.success
    assert np.allclose(res.x, _reference_blocks_solution(), atol=1e-5)


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp"])
def test_several_constraint_blocks_hess(method):
    """Check a problem with several constraint blocks and the Lagrangian Hessian."""
    res = minimize(
        _sphere,
        np.zeros(3),
        jac=_sphere_der,
        hess=_sphere_hess,
        method=method,
        tol=1e-8,
        constraints=_blocks_constraints(with_hess=True),
        options={"max_iterations": 1000},
    )
    assert res.success
    assert np.allclose(res.x, _reference_blocks_solution(), atol=1e-5)


def test_hessian_cache():
    """Check that the Hessians are evaluated once per point, not once per product."""
    points = []

    def hess(x):
        points.append(np.array(x))
        return rosen_hess(x)

    res = minimize(
        rosen,
        np.array([1.3, 0.7, 0.8]),
        jac=rosen_der,
        hess=hess,
        method="filtersqp",
        options={"max_iterations": 10000},
    )
    assert res.success
    assert np.allclose(res.x, np.ones(3), atol=1e-4)
    # Successive evaluations are at different points.
    assert all(not np.array_equal(p, q) for p, q in zip(points[:-1], points[1:]))
    # nhev counts the calls to hess, which Uno does not count for a Hessian operator.
    assert res.nhev == len(points) > 0


@pytest.mark.parametrize(
    "constraint",
    [
        NonlinearConstraint(fun=lambda x: 1.0 - 2 * x, lb=0.0, ub=np.inf),
        NonlinearConstraint(fun=lambda x: 1.0 - 2 * x, jac="2-point", lb=0.0, ub=np.inf),
        {"type": "ineq", "fun": lambda x: 1.0 - 2 * x},
        {"type": "ineq", "fun": lambda x: 1.0 - 2 * x, "jac": "2-point"},
    ],
)
def test_default_finite_difference_constraint_jacobian(constraint):
    """Check the constraint Jacobians approximated by default with finite differences."""
    res = minimize(rosen, np.zeros(2), jac=rosen_der, constraints=constraint)
    assert res.success
    assert np.allclose(res.x, [0.5, 0.25], atol=1e-4)


@pytest.mark.parametrize(
    "constraint",
    [
        NonlinearConstraint(fun=lambda x: x, jac="4-point", lb=0.0, ub=1.0),
        NonlinearConstraint(fun=lambda x: x, jac=np.eye(2), lb=0.0, ub=1.0),
        {"type": "ineq", "fun": lambda x: x, "jac": "4-point"},
        {"type": "ineq", "fun": lambda x: x, "jac": np.eye(2)},
    ],
)
def test_invalid_constraint_jacobian(constraint):
    """Check that an invalid constraint Jacobian is rejected before optimizing."""
    with pytest.raises(
        TypeError, match="constraint 0 must be a callable or a finite-difference"
    ):
        minimize(rosen, np.zeros(2), jac=rosen_der, constraints=constraint)


@pytest.mark.parametrize(
    "constraint",
    [
        NonlinearConstraint(fun=lambda x: x, jac=lambda x: np.eye(2), lb=0.0, ub=1.0),
        NonlinearConstraint(
            fun=lambda x: x, jac=lambda x: np.eye(2), hess="2-point", lb=0.0, ub=1.0
        ),
        {"type": "ineq", "fun": lambda x: x, "jac": lambda x: np.eye(2)},
    ],
)
def test_missing_constraint_hessian(constraint):
    """Check that a nonlinear constraint without Hessian is rejected when hess is given."""
    with pytest.raises(ValueError, match="constraint 1 has no callable Hessian"):
        minimize(
            rosen,
            np.zeros(2),
            jac=rosen_der,
            hess=rosen_hess,
            constraints=[LinearConstraint(np.eye(2), lb=-1.0, ub=1.0), constraint],
        )


def test_non_callable_constraint_hessian_ignored_without_hess():
    """Check that the default BFGS() Hessian of a NonlinearConstraint is ignored."""
    res = minimize(
        rosen,
        np.zeros(2),
        jac=rosen_der,
        constraints=NonlinearConstraint(
            fun=lambda x: x, jac=lambda x: np.eye(2), lb=-1.0, ub=0.5
        ),
    )
    assert res.success


def test_success_flag():
    """Check the success flag and the message when Uno does not converge."""
    res = minimize(
        rosen,
        np.zeros(2),
        jac=rosen_der,
        options={"max_iterations": 2},
    )
    assert res.success is False
    assert "ITERATION_LIMIT" in res.message

    res = minimize(rosen, np.zeros(2), jac=rosen_der)
    assert res.success is True
    assert res.message == (
        "Optimization status: SUCCESS; solution status: FEASIBLE_KKT_POINT."
    )
    assert res.status == 0
    assert type(res.status) is int
    assert res.optimization_status == unopy.OptimizationStatus.SUCCESS
    assert res.solution_status == unopy.SolutionStatus.FEASIBLE_KKT_POINT
    assert type(res.x) is np.ndarray


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp", "ipopt"])
def test_infeasible_problem(method):
    """Check that the run is not successful at an infeasible stationary point."""
    res = minimize(
        rosen,
        np.zeros(2),
        jac=rosen_der,
        method=method,
        constraints=NonlinearConstraint(
            fun=lambda x: x @ x, jac=lambda x: 2 * x, lb=-np.inf, ub=-1.0
        ),
        options={"max_iterations": 200},
    )
    assert res.success is False
    assert res.solution_status == unopy.SolutionStatus.INFEASIBLE_STATIONARY_POINT
    assert "INFEASIBLE_STATIONARY_POINT" in res.message
    # The maximum violation of x @ x <= -1, in the units of the constraint.
    assert res.maxcv == pytest.approx(res.x @ res.x + 1.0)


def _spy_set_option(monkeypatch):
    """Record the options passed to Uno."""
    options = {}
    set_option = unopy.UnoSolver.set_option

    def spy(solver, name, value):
        set_option(solver, name, value)
        options[name] = value

    monkeypatch.setattr(unopy.UnoSolver, "set_option", spy)
    return options


def test_int_tolerance_and_float_option(monkeypatch):
    """Check that an int is accepted for tol and for a float option."""
    options = _spy_set_option(monkeypatch)
    res = minimize(
        rosen,
        np.zeros(2),
        jac=rosen_der,
        tol=0,
        options={"time_limit": 100, "max_iterations": 5, "print_solution": False},
    )
    assert res.nit == 5
    assert options == {
        "primal_tolerance": 0.0,
        "dual_tolerance": 0.0,
        "time_limit": 100.0,
        "max_iterations": 5,
        "print_solution": False,
    }
    assert type(options["time_limit"]) is float
    assert type(options["max_iterations"]) is int


@pytest.mark.parametrize(
    ("name", "value", "message"),
    [
        (
            "print_solution",
            1,
            "The Uno option 'print_solution' does not accept the value 1 of type int",
        ),
        (
            "max_iterations",
            True,
            "The Uno option 'max_iterations' does not accept the value True "
            "of type bool",
        ),
        (
            "TR_radius",
            "big",
            "The Uno option 'TR_radius' does not accept the value 'big' of type str",
        ),
    ],
)
def test_wrong_option_type(name, value, message):
    """Check the error raised for a value of the wrong type.

    An int is not converted for a bool option and a bool is kept as is.
    """
    with pytest.raises(TypeError, match=message):
        minimize(rosen, np.zeros(2), jac=rosen_der, options={name: value})


@pytest.mark.parametrize("name", ["maxiter", "max_iter", "max_iterations"])
def test_max_iterations_aliases(name):
    """Check the aliases of the option max_iterations."""
    res = minimize(rosen, np.zeros(2), jac=rosen_der, options={name: 2})
    assert res.nit == 2
    assert "ITERATION_LIMIT" in res.message


@pytest.mark.parametrize("jac", ["4-point", np.zeros(2), 1])
def test_invalid_jac(jac):
    """Check that an invalid objective gradient is rejected before optimizing."""
    calls = []

    def fun(x):
        calls.append(x)
        return rosen(x)

    with pytest.raises(ValueError, match="'jac' must be a callable, True, None"):
        minimize(fun, np.zeros(2), jac=jac)
    assert not calls


def _counting(function):
    """Return a function counting its calls and the list of the calls."""
    calls = []

    def wrapped_function(x, *args):
        calls.append(np.array(x))
        return function(x, *args)

    return wrapped_function, calls


def test_constraint_size_from_bounds(monkeypatch):
    """Check that the size of a constraint is deduced from its non-scalar bounds."""
    fun, calls = _counting(lambda x: 1.0 - 2 * x)
    constraint = NonlinearConstraint(
        fun=fun, jac=lambda x: -2 * np.eye(2), lb=np.zeros(2), ub=np.inf
    )
    calls_before_optimize = []
    optimize = unopy.UnoSolver.optimize

    def spy(solver, model):
        calls_before_optimize.append(len(calls))
        return optimize(solver, model)

    monkeypatch.setattr(unopy.UnoSolver, "optimize", spy)
    res = minimize(rosen, np.zeros(2), jac=rosen_der, constraints=constraint)
    assert res.success
    assert np.allclose(res.x, [0.5, 0.25], atol=1e-4)
    # The constraint is not evaluated before the optimization.
    assert calls_before_optimize == [0]


@pytest.mark.parametrize("kind", ["NonlinearConstraint", "dict"])
def test_constraint_evaluated_once_at_x0(kind):
    """Check that a constraint is evaluated only once at x0 when its size is unknown."""
    fun, calls = _counting(lambda x: 1.0 - 2 * x)
    jac = lambda x: -2 * np.eye(2)
    if kind == "dict":
        constraint = {"type": "ineq", "fun": fun, "jac": jac}
    else:
        constraint = NonlinearConstraint(fun=fun, jac=jac, lb=0.0, ub=np.inf)
    res = minimize(rosen, np.zeros(2), jac=rosen_der, constraints=constraint)
    assert res.success
    assert np.allclose(res.x, [0.5, 0.25], atol=1e-4)
    assert sum(np.array_equal(x, np.zeros(2)) for x in calls) == 1


def test_bounds_not_mutated():
    """Check that a Bounds instance is not modified."""
    bounds = Bounds(-1.0, 0.5)
    res = minimize(rosen, np.zeros(2), jac=rosen_der, bounds=bounds)
    assert res.success
    assert np.allclose(res.x, [0.5, 0.25], atol=1e-4)
    assert bounds.lb == -1.0
    assert bounds.ub == 0.5


def test_bounds_pairs_with_none():
    """Check bounds given as a sequence of pairs with None meaning no bound."""
    bounds = [(None, 0.5), (-1, None)]
    res = minimize(rosen, np.zeros(2), jac=rosen_der, bounds=bounds)
    assert res.success
    assert np.allclose(res.x, [0.5, 0.25], atol=1e-4)
    assert bounds == [(None, 0.5), (-1, None)]


@pytest.mark.parametrize(
    ("bounds", "message"),
    [
        (Bounds(np.zeros(3), np.ones(3)), "broadcastable to the shape"),
        ([(0.0, 1.0)], "sequence of 2"),
        ([(0.0, 1.0, 2.0), (0.0, 1.0)], "sequence of 2"),
        ([(1.0, 0.0), (0.0, 1.0)], "lower bound is greater"),
        (Bounds([0.0, 2.0], [1.0, 1.0]), "lower bound is greater"),
    ],
)
def test_invalid_bounds(bounds, message):
    """Check that invalid bounds are rejected before optimizing."""
    with pytest.raises(ValueError, match=message):
        minimize(rosen, np.zeros(2), jac=rosen_der, bounds=bounds)


def test_objective_returning_array_of_size_1():
    """Check that an objective returning an array of shape (1,) raises no warning."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        res = minimize(lambda x: np.array([rosen(x)]), np.zeros(2), jac=rosen_der)
    assert res.success
    assert np.allclose(res.x, np.ones(2), atol=1e-4)


def test_objective_returning_array_of_size_2():
    """Check that an objective returning an array of size 2 is rejected."""
    with pytest.raises(
        ValueError, match=r"must return a scalar, got an array of shape \(2,\)"
    ):
        minimize(lambda x: np.array([rosen(x)] * 2), np.zeros(2), jac=rosen_der)


def test_guard_returns_zeros_after_failure():
    """Check that the guarded callbacks return zeros once a user function has failed."""
    guard = _EvaluationErrorGuard()
    calls = []

    def failing(x, output=None):
        calls.append(x)
        raise UserError("Failure.")

    objective = guard.wrap_objective(failing)
    inplace = guard.wrap_inplace(failing)
    assert objective(np.ones(2)) == 0.0
    assert guard.terminate()
    output = np.ones(3)
    inplace(np.ones(2), output)
    assert np.array_equal(output, np.zeros(3))
    # The user function is no longer called after the first failure.
    assert len(calls) == 1
    with pytest.raises(UserError, match="Failure."):
        guard.raise_exception()


def _spy_set_constraints(monkeypatch):
    """Record the number of nonzeros of the Jacobian passed to Uno.

    Returns:
        The numbers of nonzeros passed to Uno
        and the list of the calls of the Jacobian callback by Uno.
    """
    numbers_of_nonzeros = []
    uno_jacobian_calls = []
    set_constraints = unopy.Model.set_constraints

    def spy(model, *arguments):
        numbers_of_nonzeros.append(arguments[4])
        jacobian = arguments[7]

        def counted_jacobian(x, values):
            uno_jacobian_calls.append(np.array(x))
            jacobian(x, values)

        set_constraints(model, *arguments[:7], counted_jacobian)

    monkeypatch.setattr(unopy.Model, "set_constraints", spy)
    return numbers_of_nonzeros, uno_jacobian_calls


@pytest.mark.parametrize("sparse_format", ["csr", "csc", "coo"])
def test_sparse_linear_constraint(monkeypatch, sparse_format):
    """Check a LinearConstraint with a sparse matrix, kept sparse."""
    numbers_of_nonzeros, _ = _spy_set_constraints(monkeypatch)
    A = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 1.0]]))
    res = minimize(
        rosen,
        np.zeros(2),
        jac=rosen_der,
        constraints=LinearConstraint(
            A.asformat(sparse_format), ub=[0.5, 1.0, 0.5], lb=-np.inf
        ),
    )
    assert res.success
    assert np.allclose(res.x, [0.5, 0.25], atol=1e-4)
    # Only the 2 stored entries are passed to Uno, not the 3 x 2 dense matrix.
    assert numbers_of_nonzeros == [2]


def _sparse_blocks_constraints(sparse_format, kind, jacobian_calls):
    """Return the constraints of _blocks_constraints with sparse Jacobians."""

    def jac1(x):
        jacobian_calls.append(np.array(x))
        return sp.csr_matrix(_A1 + np.array([[0.0, 0.0, 2 * x[2]]])).asformat(
            sparse_format
        )

    def jac2(x):
        jacobian_calls.append(np.array(x))
        return sp.csr_matrix(_A2).asformat(sparse_format)

    def fun1(x):
        return _A1 @ x + x[2] ** 2

    def fun2(x):
        return _A2 @ x

    if kind == "dict":
        # The equivalent inequality constraints g(x) >= 0.
        return [
            {"type": "ineq", "fun": lambda x: 2.0 - fun1(x), "jac": lambda x: -jac1(x)},
            {"type": "ineq", "fun": lambda x: 1.0 - fun2(x), "jac": lambda x: -jac2(x)},
            {"type": "ineq", "fun": lambda x: 1.0 + fun2(x), "jac": jac2},
        ]
    return [
        NonlinearConstraint(fun=fun1, jac=jac1, lb=-np.inf, ub=2.0),
        NonlinearConstraint(fun=fun2, jac=jac2, lb=-1.0, ub=1.0),
    ]


@pytest.mark.parametrize("kind", ["NonlinearConstraint", "dict"])
@pytest.mark.parametrize("sparse_format", ["csr", "csc", "coo"])
def test_sparse_constraint_jacobian(monkeypatch, sparse_format, kind):
    """Check constraints with sparse Jacobians, kept sparse and not evaluated more."""
    numbers_of_nonzeros, uno_jacobian_calls = _spy_set_constraints(monkeypatch)
    jacobian_calls = []
    res = minimize(
        _sphere,
        np.zeros(3),
        jac=_sphere_der,
        method="filtersqp",
        tol=1e-8,
        constraints=_sparse_blocks_constraints(sparse_format, kind, jacobian_calls),
        options={"max_iterations": 1000},
    )
    assert res.success
    assert np.allclose(res.x, _reference_blocks_solution(), atol=1e-5)
    # The stored entries only: 3 for the block of _A1 and 5 per block of _A2,
    # instead of 3 x 3 for each block of the dense Jacobian.
    assert numbers_of_nonzeros == ([3 + 5 + 5] if kind == "dict" else [3 + 5])
    # The evaluation at x0 giving the sparsity pattern is reused by Uno.
    number_of_blocks = 3 if kind == "dict" else 2
    assert len(jacobian_calls) == number_of_blocks * len(uno_jacobian_calls)


@pytest.mark.parametrize("method", ["filtersqp", "ipopt"])
def test_dense_constraint_jacobian_calls(monkeypatch, method):
    """Check that a dense Jacobian is not evaluated more often than Uno requires."""
    _, uno_jacobian_calls = _spy_set_constraints(monkeypatch)
    jac, calls = _counting(lambda x: 2 * np.diag(x))
    res = minimize(
        rosen,
        np.array([1.3, 0.7, 0.8]),
        jac=rosen_der,
        method=method,
        tol=1e-3,
        constraints=NonlinearConstraint(
            fun=lambda x: x**2, jac=jac, lb=np.full(3, 0.1), ub=np.full(3, 0.8)
        ),
        options={"max_iterations": 10000},
    )
    assert res.success
    # The evaluation at x0 giving the sparsity pattern is reused by Uno.
    assert len(calls) == len(uno_jacobian_calls)


def test_sparse_jacobian_entry_outside_pattern():
    """Check that a nonzero entry outside the sparsity pattern at x0 is an error."""

    def jac(x):
        # The entry (0, 1) is zero at x0 and then not stored.
        jacobian = sp.csr_matrix(np.array([[1.0, 2 * x[1]]]))
        jacobian.eliminate_zeros()
        return jacobian

    with pytest.raises(
        ValueError, match=r"constraint 0 has a nonzero entry at \(0, 1\) outside"
    ):
        minimize(
            _sphere,
            np.zeros(2),
            jac=_sphere_der,
            constraints=NonlinearConstraint(
                fun=lambda x: x[0] + x[1] ** 2, jac=jac, lb=-np.inf, ub=1.0
            ),
        )


def test_jacobian_block_sparsity_pattern():
    """Check the values of a sparse Jacobian block written in its pattern order."""
    x0 = np.zeros(3)
    value_at_x0 = sp.coo_matrix(
        ([1.0, 2.0, 3.0, 4.0], ([1, 0, 1, 0], [2, 0, 0, 1])), shape=(2, 3)
    )
    values = [None]
    block = _JacobianBlock(value_at_x0, x0, 2, 0, jacobian=lambda x: values[0])
    # The stored entries in row-major order.
    assert block.nnz == 4
    assert block.rows.tolist() == [0, 0, 1, 1]
    assert block.cols.tolist() == [0, 1, 0, 2]
    output = np.empty(4)
    # The value at x0 is reused.
    block.write(x0, output)
    assert output.tolist() == [2.0, 4.0, 3.0, 1.0]

    # An entry of the pattern missing from a later evaluation is zero.
    values[0] = sp.csc_matrix(([5.0, 6.0], ([0, 1], [1, 2])), shape=(2, 3))
    block.write(np.ones(3), output)
    assert output.tolist() == [0.0, 5.0, 0.0, 6.0]

    # A dense value is accepted if its nonzeros are in the pattern.
    values[0] = np.array([[7.0, 0.0, 0.0], [8.0, 0.0, 9.0]])
    block.write(np.ones(3), output)
    assert output.tolist() == [7.0, 0.0, 8.0, 9.0]

    # An explicitly stored zero outside the pattern is ignored.
    values[0] = sp.csr_matrix(([1.0, 0.0], ([0, 0], [0, 2])), shape=(2, 3))
    block.write(np.ones(3), output)
    assert output.tolist() == [1.0, 0.0, 0.0, 0.0]

    values[0] = sp.csr_matrix(([1.0], ([1], [1])), shape=(2, 3))
    with pytest.raises(ValueError, match=r"nonzero entry at \(1, 1\) outside"):
        block.write(np.ones(3), output)

    values[0] = sp.csr_matrix((3, 3))
    with pytest.raises(ValueError, match=r"must have the shape \(2, 3\)"):
        block.write(np.ones(3), output)


# --- Hessian formats (dense, sparse, LinearOperator) ---


class _NoDenseOperator(LinearOperator):
    """A LinearOperator which must not be converted to a dense matrix."""

    def __init__(self, matrix):
        super().__init__(float, matrix.shape)
        self._matrix = sp.csr_matrix(matrix)

    def _matvec(self, vector):
        return self._matrix @ np.ravel(vector)

    def _rmatvec(self, vector):
        return self._matrix.T @ np.ravel(vector)

    def __array__(self, *args, **kwargs):
        raise AssertionError("The Hessian operator has been densified.")

    def todense(self):
        raise AssertionError("The Hessian operator has been densified.")

    def toarray(self):
        raise AssertionError("The Hessian operator has been densified.")


class _NoDenseCsr(sp.csr_matrix):
    """A sparse matrix which must not be converted to a dense matrix."""

    def __array__(self, *args, **kwargs):
        raise AssertionError("The sparse Hessian has been densified.")

    def todense(self, *args, **kwargs):
        raise AssertionError("The sparse Hessian has been densified.")

    def toarray(self, *args, **kwargs):
        raise AssertionError("The sparse Hessian has been densified.")


_HESSIAN_FORMATS = {
    "dense": lambda matrix: matrix,
    "sparse": lambda matrix: _NoDenseCsr(matrix),
    "LinearOperator": lambda matrix: _NoDenseOperator(matrix),
}


@pytest.mark.parametrize("hessian_format", sorted(_HESSIAN_FORMATS))
@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp"])
def test_hessian_formats(method, hessian_format):
    """Check that sparse and LinearOperator Hessians give the dense solution.

    The sparse matrices and the LinearOperator objects raise if densified.
    """
    to_format = _HESSIAN_FORMATS[hessian_format]
    constraints = _blocks_constraints()
    constraints[0].hess = lambda x, v: to_format(2 * v[0] * np.diag([0.0, 0.0, 1.0]))
    constraints[1].hess = lambda x, v: to_format(np.zeros((3, 3)))
    res = minimize(
        _sphere,
        np.zeros(3),
        jac=_sphere_der,
        hess=lambda x: to_format(_sphere_hess(x)),
        method=method,
        tol=1e-8,
        constraints=constraints,
        options={"max_iterations": 1000},
    )
    assert res.success
    assert np.allclose(res.x, _reference_blocks_solution(), atol=1e-5)


@pytest.mark.parametrize("hessian_format", sorted(_HESSIAN_FORMATS))
def test_hessian_formats_rosen(hessian_format):
    """Check the Rosenbrock problem with a sparse or LinearOperator Hessian."""
    to_format = _HESSIAN_FORMATS[hessian_format]
    res = minimize(
        rosen,
        np.array([1.3, 0.7, 0.8]),
        jac=rosen_der,
        hess=lambda x: to_format(rosen_hess(x)),
        options={"max_iterations": 10000},
    )
    assert res.success
    assert np.allclose(res.x, np.ones(3), atol=1e-4)


def test_lagrangian_hessian_operator_terms():
    """Check the term-by-term product and the cache of the Hessian operator."""
    x0 = np.array([1.0, 2.0])
    objective_calls = []
    constraint_calls = []

    def hess(x):
        objective_calls.append(np.array(x))
        return _NoDenseOperator(np.array([[2.0, 1.0], [1.0, 3.0]]))

    def constraint_hess(x, v):
        constraint_calls.append((np.array(x), np.array(v)))
        # The sum of v_i times the Hessian of the i-th component (x0^2, x0 * x1).
        return _NoDenseCsr(v[0] * np.diag([2.0, 0.0]) + v[1] * np.array([[0.0, 1.0], [1.0, 0.0]]))

    constraints = _normalize_constraints(
        [
            LinearConstraint(np.eye(2), lb=-1.0, ub=1.0),
            NonlinearConstraint(
                fun=lambda x: np.array([x[0] ** 2, x[0] * x[1]]),
                jac=lambda x: np.array([[2 * x[0], 0.0], [x[1], x[0]]]),
                hess=constraint_hess,
                lb=-np.inf,
                ub=1.0,
            ),
        ],
        x0,
        hess,
    )
    operator = _LagrangianHessianOperator(hess, (), constraints)

    def product(x, obj_mult, multipliers, vector):
        result = np.empty(2)
        operator(x, True, obj_mult, np.array(multipliers), np.array(vector), result)
        return result

    def expected(obj_mult, multipliers, vector):
        hessian = obj_mult * np.array([[2.0, 1.0], [1.0, 3.0]])
        hessian = hessian + multipliers[2] * np.diag([2.0, 0.0])
        hessian = hessian + multipliers[3] * np.array([[0.0, 1.0], [1.0, 0.0]])
        return hessian @ np.array(vector)

    multipliers = [5.0, 6.0, 0.5, -2.0]
    assert np.allclose(product(x0, 2.0, multipliers, [1.0, -1.0]), expected(2.0, multipliers, [1.0, -1.0]))
    # The constraint Hessian receives the multipliers of its own block.
    assert np.array_equal(constraint_calls[0][1], [0.5, -2.0])
    # Another product at the same point reuses the terms.
    assert np.allclose(product(x0, 2.0, multipliers, [0.0, 1.0]), expected(2.0, multipliers, [0.0, 1.0]))
    assert len(objective_calls) == len(constraint_calls) == 1
    # A new objective multiplier only scales the cached objective term.
    assert np.allclose(product(x0, 0.5, multipliers, [1.0, 1.0]), expected(0.5, multipliers, [1.0, 1.0]))
    assert len(objective_calls) == len(constraint_calls) == 1
    # New multipliers of the constraint: only the constraint term is recomputed.
    new_multipliers = [5.0, 6.0, 1.0, 3.0]
    assert np.allclose(product(x0, 1.0, new_multipliers, [1.0, 2.0]), expected(1.0, new_multipliers, [1.0, 2.0]))
    assert len(objective_calls) == 1
    assert len(constraint_calls) == 2
    # The multipliers of the linear constraint do not matter.
    other_multipliers = [0.0, 0.0, 1.0, 3.0]
    product(x0, 1.0, other_multipliers, [1.0, 2.0])
    assert len(objective_calls) == 1
    assert len(constraint_calls) == 2
    # A new point: both terms are recomputed.
    assert np.allclose(product(x0 + 1.0, 1.0, new_multipliers, [1.0, 2.0]), expected(1.0, new_multipliers, [1.0, 2.0]))
    assert len(objective_calls) == 2
    assert len(constraint_calls) == 3
    assert operator.evaluations == 2


# --- jac=True and finite differences ---


def test_jac_true():
    """Check that fun returning the objective and its gradient is called once per point."""
    calls = []

    def fun(x):
        calls.append(np.array(x))
        return rosen(x), rosen_der(x)

    res = minimize(fun, np.zeros(2), jac=True)
    assert res.success
    assert np.allclose(res.x, np.ones(2), atol=1e-4)
    assert res.nfev == len(calls)
    # The gradient does not call fun again at the same point.
    assert all(not np.array_equal(p, q) for p, q in zip(calls[:-1], calls[1:]))


@pytest.mark.parametrize("jac", [None, False, "2-point", "3-point", "cs"])
def test_finite_difference_objective_gradient(jac):
    """Check the objective gradient approximated with finite differences."""
    fun, calls = _counting(rosen)
    res = minimize(fun, np.zeros(2), jac=jac, tol=1e-8)
    assert res.success
    assert np.allclose(res.x, np.ones(2), atol=1e-4)
    # nfev counts the calls of the finite differences.
    assert res.nfev == len(calls)


def test_finite_difference_2_point_reuses_value():
    """Check that '2-point' reuses the objective value computed by Uno at the point."""
    fun, calls = _counting(rosen)
    res = minimize(fun, np.zeros(2), jac="2-point", options={"max_iterations": 3})
    # 2 calls per gradient for the 2 perturbed points, 1 per objective evaluation.
    assert res.nit == 3
    assert sum(np.array_equal(x, np.zeros(2)) for x in calls) == 1


def _finite_difference_blocks_constraints(method, kind):
    """Return the constraints of _blocks_constraints with finite-difference Jacobians."""

    def fun1(x):
        return _A1 @ x + x[2] ** 2

    def fun2(x):
        return _A2 @ x

    if kind == "dict":
        jac = {} if method is None else {"jac": method}
        return [
            {"type": "ineq", "fun": lambda x: 2.0 - fun1(x), **jac},
            {"type": "ineq", "fun": lambda x: 1.0 - fun2(x), **jac},
            {"type": "ineq", "fun": lambda x: 1.0 + fun2(x), **jac},
        ]
    jac = {} if method is None else {"jac": method}
    return [
        NonlinearConstraint(fun=fun1, lb=-np.inf, ub=2.0, **jac),
        NonlinearConstraint(fun=fun2, lb=-1.0, ub=1.0, **jac),
    ]


@pytest.mark.parametrize("kind", ["NonlinearConstraint", "dict"])
@pytest.mark.parametrize("method", [None, "2-point", "3-point", "cs"])
def test_finite_difference_constraint_jacobian(monkeypatch, method, kind):
    """Check constraint Jacobians approximated with finite differences (dense)."""
    numbers_of_nonzeros, _ = _spy_set_constraints(monkeypatch)
    res = minimize(
        _sphere,
        np.zeros(3),
        jac=_sphere_der,
        tol=1e-8,
        constraints=_finite_difference_blocks_constraints(method, kind),
        options={"max_iterations": 1000},
    )
    assert res.success
    assert np.allclose(res.x, _reference_blocks_solution(), atol=1e-5)
    # The sparsity pattern of a finite-difference Jacobian is dense.
    assert numbers_of_nonzeros == ([5 * 3] if kind == "dict" else [3 * 3])


@pytest.mark.parametrize("method", ["filtersqp", "ipopt"])
@pytest.mark.parametrize("failing", ["fun", "constraint"])
def test_finite_difference_exception(method, failing):
    """Check that an exception raised during the finite differences is re-raised."""
    functions = {"fun": rosen, "constraint": lambda x: x**2}
    # The 4th call is a perturbed point of the first finite differences.
    functions[failing], calls = raise_after(functions[failing], 3)
    with pytest.raises(UserError, match="Failure at call 4."):
        minimize(
            functions["fun"],
            np.array([1.3, 0.7, 0.8]),
            jac="2-point",
            method=method,
            constraints=NonlinearConstraint(
                fun=functions["constraint"], lb=np.full(3, 0.1), ub=np.full(3, 0.8)
            ),
            options={"max_iterations": 10000},
        )
    assert len(calls) == 4


@pytest.mark.parametrize(
    ("hess", "hessian_model"), [(BFGS(), "LBFGS"), (SR1(), "LSR1")]
)
@pytest.mark.parametrize("constrained", [False, True])
def test_quasi_newton_hessian(monkeypatch, hess, hessian_model, constrained):
    """Check that BFGS() and SR1() select the quasi-Newton models of Uno."""
    options = _spy_set_option(monkeypatch)
    set_operator = []
    monkeypatch.setattr(
        unopy.Model,
        "set_lagrangian_hessian_operator",
        lambda *args: set_operator.append(args),
    )
    if constrained:
        # A callable constraint Hessian is ignored with a quasi-Newton model.
        constraints = _blocks_constraints(with_hess=True)
        constraints[0].hess = BFGS()
        reference = _reference_blocks_solution()
        res = minimize(
            _sphere, np.zeros(3), jac=_sphere_der, hess=hess, tol=1e-8,
            constraints=constraints, options={"max_iterations": 1000},
        )
    else:
        reference = np.ones(2)
        res = minimize(
            rosen, np.zeros(2), jac=rosen_der, hess=hess,
            options={"max_iterations": 1000},
        )
    assert res.success
    assert np.allclose(res.x, reference, atol=1e-5)
    assert options["hessian_model"] == hessian_model
    assert not set_operator


def test_quasi_newton_hessian_overridden_by_option(monkeypatch):
    """Check that the option hessian_model overrides the one set from hess."""
    options = []
    set_option = unopy.UnoSolver.set_option

    def spy(solver, name, value):
        set_option(solver, name, value)
        options.append((name, value))

    monkeypatch.setattr(unopy.UnoSolver, "set_option", spy)
    res = minimize(
        rosen, np.zeros(2), jac=rosen_der, hess=SR1(), options={"hessian_model": "LBFGS"}
    )
    assert res.success
    assert options == [("hessian_model", "LSR1"), ("hessian_model", "LBFGS")]


def test_quasi_newton_hessian_filterslp(monkeypatch):
    """Check that the Hessian model is not set for the SLP method."""
    options = _spy_set_option(monkeypatch)
    res = minimize(
        rosen, np.zeros(2), jac=rosen_der, hess=BFGS(), method="filterslp",
        options={"max_iterations": 5},
    )
    assert res.nit == 5
    assert "hessian_model" not in options


@pytest.mark.parametrize(
    ("hess", "error", "message"),
    [
        ("2-point", NotImplementedError, "finite-difference Hessian '2-point'"),
        ("cs", NotImplementedError, "finite-difference Hessian 'cs'"),
        (np.eye(2), ValueError, "'hess' must be a callable, None, BFGS"),
        ("exact", ValueError, "'hess' must be a callable, None, BFGS"),
    ],
)
def test_invalid_hess(hess, error, message):
    """Check that an unsupported Hessian is rejected before optimizing."""
    fun, calls = _counting(rosen)
    with pytest.raises(error, match=message):
        minimize(fun, np.zeros(2), jac=rosen_der, hess=hess)
    assert not calls


def test_hessian_update_strategy_constraint_with_callable_hess():
    """Check that a constraint Hessian update strategy is rejected with a callable hess."""
    with pytest.raises(ValueError, match="constraint 0 has no callable Hessian"):
        minimize(
            rosen,
            np.zeros(2),
            jac=rosen_der,
            hess=rosen_hess,
            constraints=NonlinearConstraint(
                fun=lambda x: x, jac=lambda x: np.eye(2), hess=SR1(), lb=0.0, ub=1.0
            ),
        )


# --- callback ---


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp", "ipopt"])
def test_callback_intermediate_result(method):
    """Check the callback called with an OptimizeResult at each accepted iterate."""
    results = []

    def callback(intermediate_result):
        results.append(intermediate_result)

    fun, calls = _counting(_sphere)
    constraints = _blocks_constraints()
    res = minimize(
        fun, np.zeros(3), jac=_sphere_der, method=method, tol=1e-8,
        constraints=constraints, bounds=Bounds(-10.0, 10.0), callback=callback,
        options={"max_iterations": 1000},
    )
    assert res.success
    # One call per accepted iterate.
    assert len(results) == res.nit
    assert [result.nit for result in results] == list(range(1, res.nit + 1))
    for result in results:
        assert isinstance(result, OptimizeResult)
        assert result.x.shape == (3,)
        assert result.fun == pytest.approx(_sphere(result.x))
        assert result.constraint_multipliers.shape == (3,)
        assert result.lower_bound_multipliers.shape == (3,)
        assert result.upper_bound_multipliers.shape == (3,)
        assert isinstance(result.stationarity, float)
        assert isinstance(result.primal_feasibility, float)
        assert isinstance(result.complementarity, float)
        assert isinstance(result.objective_multiplier, float)
    assert np.allclose(results[-1].x, res.x)
    assert results[-1].fun == pytest.approx(res.fun)
    # The multipliers at the last iterate are those of the result.
    assert np.allclose(results[-1].constraint_multipliers, np.concatenate(res.v[:2]), atol=1e-6)
    # The callback makes no additional call to fun.
    fun_without_callback, calls_without_callback = _counting(_sphere)
    minimize(
        fun_without_callback, np.zeros(3), jac=_sphere_der, method=method, tol=1e-8,
        constraints=constraints, bounds=Bounds(-10.0, 10.0),
        options={"max_iterations": 1000},
    )
    assert len(calls) == len(calls_without_callback)


def test_callback_legacy_signature():
    """Check the callback called with a copy of the iterate."""
    iterates = []
    res = minimize(rosen, np.zeros(2), jac=rosen_der, callback=iterates.append)
    assert res.success
    assert len(iterates) == res.nit
    assert all(isinstance(x, np.ndarray) and x.shape == (2,) for x in iterates)
    assert all(not np.shares_memory(x, y) for x, y in zip(iterates[:-1], iterates[1:]))
    assert np.allclose(iterates[-1], res.x)


@pytest.mark.parametrize("method", ["filtersqp", "ipopt"])
def test_callback_stop_iteration(method):
    """Check that StopIteration raised by the callback stops Uno without error."""
    iterates = []

    def callback(xk):
        iterates.append(xk)
        if len(iterates) == 3:
            raise StopIteration

    res = minimize(rosen, np.zeros(2), jac=rosen_der, method=method, callback=callback)
    assert len(iterates) == 3
    assert res.nit == 3
    assert res.success is False
    assert res.message == "`callback` raised `StopIteration`."
    assert res.status == 99
    assert np.allclose(res.x, iterates[-1])


def test_callback_exception():
    """Check that an exception raised by the callback stops Uno and is re-raised."""
    iterates = []
    fun, calls = _counting(rosen)
    number_of_calls = []

    def callback(intermediate_result):
        iterates.append(intermediate_result.x)
        if len(iterates) == 2:
            number_of_calls.append(len(calls))
            raise UserError("Callback failure.")

    with pytest.raises(UserError, match="Callback failure."):
        minimize(fun, np.zeros(2), jac=rosen_der, callback=callback)
    assert len(iterates) == 2
    # Uno is stopped and fun is no longer called.
    assert number_of_calls == [len(calls)]


def test_callback_not_called_after_user_exception():
    """Check that the callback is not called once a user function has failed."""
    iterates = []
    fun, calls = raise_after(rosen, 3)
    with pytest.raises(UserError, match="Failure at call 4."):
        minimize(fun, np.zeros(2), jac=rosen_der, callback=iterates.append)
    assert len(iterates) <= 3


# --- multipliers ---


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp", "ipopt"])
@pytest.mark.parametrize("kind", ["NonlinearConstraint", "LinearConstraint", "dict"])
def test_inequality_multiplier(method, kind):
    """Check the multiplier of min x0^2 + x1^2 s.t. x0 + x1 >= 1.

    The KKT conditions 2 x + v (1, 1) = 0 at x = (0.5, 0.5) give v = -1.
    """
    if kind == "NonlinearConstraint":
        constraint = NonlinearConstraint(
            fun=lambda x: x[0] + x[1], jac=lambda x: np.ones((1, 2)), lb=1.0, ub=np.inf
        )
    elif kind == "LinearConstraint":
        constraint = LinearConstraint(np.ones((1, 2)), lb=1.0, ub=np.inf)
    else:
        constraint = {
            "type": "ineq", "fun": lambda x: x[0] + x[1] - 1.0, "jac": lambda x: np.ones(2)
        }
    res = minimize(
        lambda x: x @ x, np.zeros(2), jac=lambda x: 2 * x, method=method,
        constraints=constraint, tol=1e-8,
    )
    assert res.success
    assert np.allclose(res.x, [0.5, 0.5], atol=1e-6)
    assert len(res.v) == 1
    assert np.allclose(res.v[0], [-1.0], atol=1e-6)
    assert np.allclose(res.constr[0], [1.0 if kind != "dict" else 0.0], atol=1e-6)
    assert np.allclose(res.lower_bound_multipliers, 0.0, atol=1e-6)
    assert np.allclose(res.upper_bound_multipliers, 0.0, atol=1e-6)


@pytest.mark.parametrize("method", ["filtersqp", "ipopt"])
def test_upper_bound_constraint_multiplier(method):
    """Check the positive multiplier of an active upper bound of a constraint."""
    # min (x0 - 2)^2 + (x1 - 2)^2 s.t. x0 + x1 <= 2: x = (1, 1), v = 2.
    res = minimize(
        lambda x: ((x - 2.0) ** 2).sum(), np.zeros(2), jac=lambda x: 2 * (x - 2.0),
        method=method, constraints=LinearConstraint(np.ones((1, 2)), ub=2.0), tol=1e-8,
    )
    assert res.success
    assert np.allclose(res.x, [1.0, 1.0], atol=1e-6)
    assert np.allclose(res.v[0], [2.0], atol=1e-6)


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp", "ipopt"])
def test_equality_and_bound_multipliers(method):
    """Check the multipliers of an equality constraint and of active bounds.

    min (x0 - 2)^2 + (x1 + 2)^2 + x2^2 + x3^2
    s.t. x2 + 2 x3 = 1, -1 <= x0, x1 <= 1;
    x = (1, -1, 0.2, 0.4), v = -0.4 for the equality,
    0 = 2 (x0 - 2) + z0 gives z0 = 2 for the upper bound of x0,
    0 = 2 (x1 + 2) + z1 gives z1 = -2 for the lower bound of x1.
    """
    res = minimize(
        lambda x: (x[0] - 2) ** 2 + (x[1] + 2) ** 2 + x[2] ** 2 + x[3] ** 2,
        np.zeros(4),
        jac=lambda x: np.array([2 * (x[0] - 2), 2 * (x[1] + 2), 2 * x[2], 2 * x[3]]),
        method=method,
        bounds=Bounds([-1.0, -1.0, -np.inf, -np.inf], [1.0, 1.0, np.inf, np.inf]),
        constraints=[
            {"type": "eq", "fun": lambda x: x[2] + 2 * x[3] - 1.0,
             "jac": lambda x: np.array([0.0, 0.0, 1.0, 2.0])},
        ],
        tol=1e-8,
    )
    assert res.success
    assert np.allclose(res.x, [1.0, -1.0, 0.2, 0.4], atol=1e-6)
    # The multipliers of the constraint, then those of the bounds.
    assert len(res.v) == 2
    assert np.allclose(res.v[0], [-0.4], atol=1e-6)
    assert np.allclose(res.v[1], [2.0, -2.0, 0.0, 0.0], atol=1e-6)
    assert np.allclose(res.upper_bound_multipliers, [2.0, 0.0, 0.0, 0.0], atol=1e-6)
    assert np.allclose(res.lower_bound_multipliers, [0.0, -2.0, 0.0, 0.0], atol=1e-6)
    # The gradient of the Lagrangian is zero.
    gradient = np.array([2 * (res.x[0] - 2), 2 * (res.x[1] + 2), 2 * res.x[2], 2 * res.x[3]])
    assert np.allclose(gradient + res.v[0] * np.array([0.0, 0.0, 1.0, 2.0]) + res.v[1], 0.0, atol=1e-6)


def test_multipliers_of_several_constraints():
    """Check that the multipliers are split per constraint in the user's order."""
    res = minimize(
        _sphere, np.zeros(3), jac=_sphere_der, tol=1e-8, constraints=_blocks_constraints(),
    )
    reference = scipy_minimize(
        _sphere, np.zeros(3), jac=_sphere_der, hess=lambda x: 2 * np.eye(3),
        method="trust-constr", constraints=_blocks_constraints(with_hess=True),
        options={"gtol": 1e-10, "xtol": 1e-12},
    )
    assert res.success
    assert [v.shape for v in res.v] == [(1,), (2,)]
    for v, v_reference in zip(res.v, reference.v):
        assert np.allclose(v, v_reference, atol=1e-3)
    assert [c.shape for c in res.constr] == [(1,), (2,)]


# --- filterslp ---


def test_filterslp_linear_problem():
    """Check that the SLP method solves a linear problem at a vertex."""
    res = minimize(
        lambda x: -x[0] - 2 * x[1],
        np.zeros(2),
        jac=lambda x: np.array([-1.0, -2.0]),
        method="filterslp",
        bounds=Bounds([0.0, 0.0], [3.0, 3.0]),
        constraints=LinearConstraint(
            [[1.0, 1.0], [1.0, -1.0]], lb=[-np.inf, -2.0], ub=[4.0, np.inf]
        ),
    )
    assert res.success
    assert np.allclose(res.x, [1.0, 3.0], atol=1e-8)
    assert res.fun == pytest.approx(-7.0)


@pytest.mark.parametrize("kind", ["NonlinearConstraint", "dict"])
def test_filterslp_nonlinear_constraint(kind):
    """Check that the SLP method solves a linear objective over a disk."""
    if kind == "dict":
        constraint = {
            "type": "ineq",
            "fun": lambda x: 1.0 - x @ x,
            "jac": lambda x: -2 * x.reshape(1, -1),
        }
    else:
        constraint = NonlinearConstraint(
            fun=lambda x: x @ x, lb=-np.inf, ub=1.0, jac=lambda x: 2 * x.reshape(1, -1)
        )
    res = minimize(
        lambda x: x[0] + x[1],
        np.array([0.5, 0.0]),
        jac=lambda x: np.ones(2),
        method="filterslp",
        tol=1e-8,
        constraints=constraint,
        options={"max_iterations": 1000},
    )
    assert res.success
    assert np.allclose(res.x, [-(0.5**0.5), -(0.5**0.5)], atol=1e-5)


# --- invalid inputs ---


class _CustomHessianUpdate(HessianUpdateStrategy):
    """A Hessian update strategy which is neither BFGS nor SR1."""


@pytest.mark.parametrize(
    ("kwargs", "error", "message"),
    [
        ({"method": "slsqp"}, ValueError, r"Unknown method 'slsqp'\. Must be one of"),
        (
            {"hess": _CustomHessianUpdate()},
            NotImplementedError,
            r"The Hessian update strategy .* is not supported; use BFGS\(\) or SR1\(\)",
        ),
        (
            {"constraints": {"type": "lt", "fun": lambda x: x, "jac": lambda x: np.eye(2)}},
            ValueError,
            "Unknown constraint type 'lt'. Must be 'eq' or 'ineq'",
        ),
        (
            {"constraints": [LinearConstraint(np.eye(2), lb=-1.0, ub=1.0), (0.0, 1.0)]},
            TypeError,
            "Unsupported constraint type: <class 'tuple'>",
        ),
        (
            {
                "constraints": NonlinearConstraint(
                    fun=lambda x: x[0], lb=-np.inf, ub=1.0, jac=lambda x: np.ones(3)
                )
            },
            ValueError,
            r"constraint 0 must have the shape \(1, 2\), got an array of shape \(3,\)",
        ),
    ],
    ids=["method", "hess", "dict_type", "constraint_type", "jacobian_size"],
)
def test_invalid_input(kwargs, error, message):
    """Check that an invalid input is rejected before optimizing."""
    fun, calls = _counting(rosen)
    with pytest.raises(error, match=message):
        minimize(fun, np.zeros(2), jac=rosen_der, **kwargs)
    assert not calls


# --- non-canonical sparse Jacobians ---


def _non_canonical_jacobian(x, sparse_format):
    """Return the Jacobian of _non_canonical_constraint with duplicate entries.

    The Jacobian [[1, 4 x1, 0], [0, 0, 3]] is stored with the entries of the row 0
    in unsorted column order, the entry (0, 1) split into two duplicates
    and the entry (1, 2) split into two duplicates.
    """
    data = np.array([2 * x[1], 1.0, 2 * x[1], 1.0, 2.0])
    rows = np.array([0, 0, 0, 1, 1])
    cols = np.array([1, 0, 1, 2, 2])
    if sparse_format == "coo":
        return sp.coo_matrix((data, (rows, cols)), shape=(2, 3))
    jacobian = sp.csr_matrix((data, cols, np.array([0, 3, 5])), shape=(2, 3))
    assert not jacobian.has_canonical_format
    return jacobian


def _non_canonical_constraint(x):
    return np.array([x[0] + 2 * x[1] ** 2, 3 * x[2]])


@pytest.mark.parametrize("sparse_format", ["csr", "coo"])
def test_jacobian_block_non_canonical(sparse_format):
    """Check a sparse Jacobian block with unsorted indices and duplicate entries."""
    x0 = np.zeros(3)
    jacobians = []

    def jac(x):
        jacobians.append(_non_canonical_jacobian(x, sparse_format))
        return jacobians[-1]

    block = _JacobianBlock(jac(x0), x0, 2, 0, jacobian=jac)
    # The duplicates are summed, including the explicitly stored zeros at x0.
    assert block.nnz == 3
    assert block.rows.tolist() == [0, 0, 1]
    assert block.cols.tolist() == [0, 1, 2]
    output = np.empty(3)
    block.write(x0, output)
    assert output.tolist() == [1.0, 0.0, 3.0]
    block.write(np.array([0.0, 1.5, 0.0]), output)
    assert output.tolist() == [1.0, 6.0, 3.0]
    # The matrices of the user are not modified.
    assert all(jacobian.nnz == 5 for jacobian in jacobians)


@pytest.mark.parametrize("sparse_format", ["csr", "coo"])
def test_non_canonical_sparse_constraint_jacobian(sparse_format):
    """Check the solution with a sparse Jacobian with unsorted and duplicate entries."""

    def solve(jac):
        return minimize(
            _sphere,
            np.zeros(3),
            jac=_sphere_der,
            tol=1e-8,
            constraints=NonlinearConstraint(
                fun=_non_canonical_constraint, lb=-np.inf, ub=[2.0, 1.0], jac=jac
            ),
        )

    res = solve(lambda x: _non_canonical_jacobian(x, sparse_format))
    reference = solve(lambda x: _non_canonical_jacobian(x, sparse_format).toarray())
    assert res.success
    assert res.nit > 1
    assert np.allclose(res.x, reference.x, atol=1e-8)
    assert np.all(_non_canonical_constraint(res.x) <= [2.0 + 1e-8, 1.0 + 1e-8])


def test_sparse_jacobian_of_dense_block():
    """Check a sparse Jacobian returned for a block whose Jacobian at x0 is dense."""
    x0 = np.zeros(3)

    def jac(x):
        jacobian = _non_canonical_jacobian(x, "csr")
        return jacobian.toarray() if np.array_equal(x, x0) else jacobian

    block = _JacobianBlock(jac(x0), x0, 2, 0, jacobian=jac)
    # The dense pattern.
    assert block.nnz == 6
    output = np.empty(6)
    block.write(np.array([0.0, 1.5, 0.0]), output)
    assert output.tolist() == [1.0, 6.0, 0.0, 0.0, 0.0, 3.0]

    def solve(jac):
        return minimize(
            _sphere,
            x0,
            jac=_sphere_der,
            tol=1e-8,
            constraints=NonlinearConstraint(
                fun=_non_canonical_constraint, lb=-np.inf, ub=[2.0, 1.0], jac=jac
            ),
        )

    res = solve(jac)
    reference = solve(lambda x: _non_canonical_jacobian(x, "csr").toarray())
    assert res.success
    assert np.allclose(res.x, reference.x, atol=1e-8)


# --- objective, finite differences and Hessian operator ---


def test_objective_jac_true_gradient_at_another_point():
    """Check jac=True when the gradient is requested at another point than the last one."""
    fun, calls = _counting(lambda x: (rosen(x), rosen_der(x)))
    objective = _Objective(fun, True, ())
    x1 = np.array([0.5, 0.5])
    x2 = np.array([-1.0, 2.0])
    assert objective.evaluate(x1) == rosen(x1)
    gradient = np.empty(2)
    # fun is called at the new point.
    objective.compute_gradient(x2, gradient)
    assert np.array_equal(gradient, rosen_der(x2))
    assert len(calls) == 2
    # The value at this point is then reused.
    assert objective.evaluate(x2) == rosen(x2)
    assert len(calls) == 2
    # Another point requires a new call.
    assert objective.evaluate(x1) == rosen(x1)
    assert len(calls) == 3
    assert objective.evaluations == 3


def test_approximate_jacobian_2_point_without_value():
    """Check the '2-point' scheme evaluating the function at the point itself."""
    function, calls = _counting(lambda x: np.array([x[0] ** 2, x[0] * x[1]]))
    x = np.array([1.0, -2.0])
    jacobian = _approximate_jacobian(function, x, "2-point")
    assert np.allclose(jacobian, [[2.0, 0.0], [-2.0, 1.0]], atol=1e-6)
    # The point itself and the 2 perturbed points.
    assert len(calls) == 3
    assert np.array_equal(calls[0], x)


def test_lagrangian_hessian_operator_same_point():
    """Check that the Hessian term is reused when Uno has not changed the point."""
    calls = []

    def hess(x):
        calls.append(np.array(x))
        return np.diag(x)

    operator = _LagrangianHessianOperator(hess, (), [])
    result = np.empty(2)
    operator(np.array([1.0, 2.0]), True, 1.0, np.empty(0), np.ones(2), result)
    assert result.tolist() == [1.0, 2.0]
    operator(np.array([1.0, 2.0]), False, 2.0, np.empty(0), np.ones(2), result)
    assert result.tolist() == [2.0, 4.0]
    assert len(calls) == 1


def test_lagrangian_hessian_operator_same_point_new_multipliers():
    """Check that a constraint term is recomputed when only the multipliers change.

    Uno may keep the point (``evaluate_at_x=False``) and change the multipliers:
    the objective term is reused but not the constraint term.
    """
    x0 = np.array([1.0, 2.0])
    objective_calls = []
    constraint_calls = []

    def hess(x):
        objective_calls.append(np.array(x))
        return np.zeros((2, 2))

    def constraint_hess(x, v):
        constraint_calls.append(np.array(v))
        return v[0] * np.diag([2.0, 2.0])

    constraints = _normalize_constraints(
        NonlinearConstraint(
            fun=lambda x: np.array([x @ x]),
            jac=lambda x: np.atleast_2d(2 * x),
            hess=constraint_hess,
            lb=-np.inf,
            ub=1.0,
        ),
        x0,
        hess,
    )
    operator = _LagrangianHessianOperator(hess, (), constraints)
    result = np.empty(2)
    operator(x0, True, 1.0, np.array([1.0]), np.ones(2), result)
    assert result.tolist() == [2.0, 2.0]
    # Same point and multipliers: both terms are reused.
    operator(x0, False, 1.0, np.array([1.0]), np.ones(2), result)
    assert result.tolist() == [2.0, 2.0]
    assert len(objective_calls) == len(constraint_calls) == 1
    # Same point, new multipliers: only the constraint term is recomputed.
    operator(x0, False, 1.0, np.array([3.0]), np.ones(2), result)
    assert result.tolist() == [6.0, 6.0]
    assert len(objective_calls) == 1
    assert len(constraint_calls) == 2


# --- callback with a signature that cannot be inspected ---


class _UninspectableCallback:
    """A callback whose signature cannot be inspected."""

    def __init__(self):
        self.calls = []

    @property
    def __signature__(self):
        raise ValueError("No signature.")

    def __call__(self, *args):
        self.calls.append(args)


def test_callback_uninspectable_signature():
    """Check that a callback whose signature cannot be inspected receives the iterate."""
    callback = _UninspectableCallback()
    with pytest.raises(ValueError, match="No signature."):
        inspect.signature(callback)
    res = minimize(rosen, np.zeros(2), jac=rosen_der, callback=callback)
    assert res.success
    assert len(callback.calls) == res.nit
    # The legacy form callback(xk).
    assert all(
        len(args) == 1 and isinstance(args[0], np.ndarray) and args[0].shape == (2,)
        for args in callback.calls
    )
    assert np.allclose(callback.calls[-1][0], res.x)


def test_hess_with_ipopt():
    """Check that a callable Hessian with 'ipopt' raises before the optimization."""
    with pytest.raises(ValueError, match="requires an explicit Hessian matrix"):
        minimize(rosen, np.zeros(2), jac=rosen_der, hess=rosen_hess, method="ipopt")


@pytest.mark.parametrize(
    ("method", "options"),
    [
        ("filterslp", {}),
        ("filtersqp", {"hessian_model": "LBFGS"}),
        ("ipopt", {"hessian_model": "LBFGS"}),
    ],
)
def test_hess_not_used(method, options):
    """Check that a callable Hessian unused by the method is ignored.

    The constraint Hessians are then not required.
    """
    calls = []

    def hess(x):
        calls.append(x)
        return rosen_hess(x)

    res = minimize(
        rosen,
        np.array([0.5, 0.5]),
        jac=rosen_der,
        hess=hess,
        method=method,
        bounds=[(-2.0, 2.0)] * 2,
        constraints={"type": "ineq", "fun": lambda x: 1.0 - x @ x},
        options={"max_iterations": 500, **options},
    )
    assert not calls
    assert res.nhev == 0


def test_lagrangian_hessian_operator_zero_multipliers():
    """Check that the terms with zero multipliers are not computed."""
    x0 = np.array([1.0, 2.0])
    objective_calls = []
    constraint_calls = []

    def hess(x):
        objective_calls.append(np.array(x))
        return np.eye(2)

    def constraint_hess(x, v):
        constraint_calls.append(np.array(v))
        return v[0] * np.eye(2)

    constraints = _normalize_constraints(
        NonlinearConstraint(
            fun=lambda x: np.array([x @ x]),
            jac=lambda x: np.atleast_2d(2 * x),
            hess=constraint_hess,
            lb=-np.inf,
            ub=1.0,
        ),
        x0,
        hess,
    )
    operator = _LagrangianHessianOperator(hess, (), constraints)
    result = np.empty(2)
    # E.g. the feasibility restoration phase: the objective multiplier is zero.
    operator(x0, True, 0.0, np.array([2.0]), np.ones(2), result)
    assert result.tolist() == [2.0, 2.0]
    assert not objective_calls
    assert len(constraint_calls) == 1
    # An inactive constraint: its multipliers are zero.
    operator(x0, True, 1.0, np.array([0.0]), np.ones(2), result)
    assert result.tolist() == [1.0, 1.0]
    assert len(objective_calls) == 1
    assert len(constraint_calls) == 1
    assert operator.evaluations == 1


def test_dense_linear_constraint_pattern(monkeypatch):
    """Check that the sparsity pattern of a dense LinearConstraint is its nonzeros."""
    numbers_of_nonzeros, _ = _spy_set_constraints(monkeypatch)
    res = minimize(
        rosen,
        np.zeros(3),
        jac=rosen_der,
        constraints=LinearConstraint(
            np.array([[1.0, 0.0, 1.0], [0.0, 0.0, 0.0], [0.0, 2.0, 0.0]]), -1.0, 0.5
        ),
    )
    assert numbers_of_nonzeros == [3]
    assert res.success
    reference = scipy_minimize(
        rosen,
        np.zeros(3),
        jac=rosen_der,
        constraints=LinearConstraint(
            np.array([[1.0, 0.0, 1.0], [0.0, 0.0, 0.0], [0.0, 2.0, 0.0]]), -1.0, 0.5
        ),
        method="trust-constr",
    )
    assert np.allclose(res.x, reference.x, atol=1e-4)


def test_jacobian_block_same_structure():
    """Check the values of a sparse Jacobian with the structure of the pattern.

    The values are copied directly; a Jacobian with another structure,
    e.g. without an explicit zero of the pattern, is still written correctly.
    """
    x0 = np.zeros(2)
    pattern = sp.csr_matrix(
        (np.array([1.0, 0.0, 3.0]), np.array([0, 1, 1]), np.array([0, 2, 3])),
        shape=(2, 2),
    )
    values = [pattern]
    block = _JacobianBlock(pattern, x0, 2, 0, jacobian=lambda x: values[0])
    assert block.nnz == 3
    output = np.empty(3)
    block.write(x0, output)
    assert output.tolist() == [1.0, 0.0, 3.0]
    values[0] = sp.csr_matrix(
        (np.array([4.0, 5.0, 6.0]), np.array([0, 1, 1]), np.array([0, 2, 3])),
        shape=(2, 2),
    )
    block.write(x0, output)
    assert output.tolist() == [4.0, 5.0, 6.0]
    values[0] = sp.csr_matrix(np.array([[7.0, 0.0], [0.0, 8.0]]))
    block.write(x0, output)
    assert output.tolist() == [7.0, 0.0, 8.0]


def test_njev_and_maxcv_with_bounds():
    """Check njev and the maximum violation of the bounds."""
    gradient_calls = []

    def jac(x):
        gradient_calls.append(x)
        return rosen_der(x)

    res = minimize(rosen, np.zeros(2), jac=jac, bounds=[(-1.0, 0.5), (-1.0, 0.5)])
    assert res.njev == len(gradient_calls)
    assert res.maxcv == pytest.approx(
        max(0.0, np.max(res.x - 0.5), np.max(-1.0 - res.x))
    )
    assert res.maxcv < 1e-6


def test_documented_options_are_uno_options():
    """Check that the options documented by minimize are the options of Uno.

    The options of Uno are read from uno/options/Options.cpp,
    except "preset", which is set from the method.
    """
    options_file = (
        Path(__file__).resolve().parents[3] / "uno" / "options" / "Options.cpp"
    )
    if not options_file.is_file():
        pytest.skip("The source of Uno is not available.")
    types = {"DOUBLE": "float", "INTEGER": "int", "BOOL": "bool", "STRING": "str"}
    uno_options = {
        name: types[kind]
        for name, kind in re.findall(
            r'\{"(\w+)",\s*OptionType::(\w+)\}', options_file.read_text()
        )
    }
    documented_options = dict(
        re.findall(r'^\s*"(\w+)"[^:\n]*: (\w+)\s*$', minimize.__doc__, re.MULTILINE)
    )
    del uno_options["preset"]
    assert documented_options == uno_options
