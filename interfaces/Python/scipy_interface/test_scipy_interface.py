# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

import numpy as np
import pytest
from scipy.optimize import rosen, rosen_der, rosen_hess, NonlinearConstraint, LinearConstraint

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
        constr = LinearConstraint(A=np.eye(2), ub=0.5 * np.zeros(2), lb=-float("Inf"))
    else:
        constr = {
            "type": "ineq",
            "fun": lambda x: 1.0 - 2 * x,
            "jac": lambda x: -2 * np.eye(2),
        }
    res = minimize(rosen, jac=rosen_der, x0=np.zeros(2), constraints=(constr,))

    assert res.success

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
