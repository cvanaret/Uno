# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

import numpy as np
import pytest
import unopy
from scipy.optimize import rosen, rosen_der, rosen_hess, Bounds, NonlinearConstraint, LinearConstraint
from scipy.optimize import minimize as scipy_minimize

from scipy_interface.scipy_uno import minimize


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


@pytest.mark.parametrize(
    "constraint",
    [
        NonlinearConstraint(fun=lambda x: x, lb=0.0, ub=1.0),
        NonlinearConstraint(fun=lambda x: x, jac="2-point", lb=0.0, ub=1.0),
        {"type": "ineq", "fun": lambda x: x},
        {"type": "ineq", "fun": lambda x: x, "jac": "2-point"},
    ],
)
def test_non_callable_constraint_jacobian(constraint):
    """Check that a non-callable constraint Jacobian is rejected before optimizing."""
    with pytest.raises(TypeError, match="constraint 0 must be a callable"):
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
    assert res.message == unopy.OptimizationStatus.ITERATION_LIMIT

    res = minimize(rosen, np.zeros(2), jac=rosen_der)
    assert res.success is True
    assert res.message == unopy.OptimizationStatus.SUCCESS


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
        ("print_solution", 1, "print_solution is not of type int"),
        ("max_iterations", True, "max_iterations is not of type bool"),
    ],
)
def test_wrong_option_type(name, value, message):
    """Check that an int is not converted for a bool option and a bool is kept as is."""
    with pytest.raises(TypeError, match=message):
        minimize(rosen, np.zeros(2), jac=rosen_der, options={name: value})


@pytest.mark.parametrize("name", ["maxiter", "max_iter", "max_iterations"])
def test_max_iterations_aliases(name):
    """Check the aliases of the option max_iterations."""
    res = minimize(rosen, np.zeros(2), jac=rosen_der, options={name: 2})
    assert res.nit == 2
    assert res.message == unopy.OptimizationStatus.ITERATION_LIMIT


@pytest.mark.parametrize("jac", [True, "2-point"])
def test_non_callable_jac(jac):
    """Check that a non-callable objective gradient is rejected before optimizing."""
    calls = []

    def fun(x):
        calls.append(x)
        return rosen(x), rosen_der(x)

    with pytest.raises(ValueError, match="'jac' must be a callable"):
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
