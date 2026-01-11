import numpy as np
import pytest
from scipy.optimize import rosen, rosen_der, NonlinearConstraint, LinearConstraint

from scipy_interface import minimize


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp", "filterslp"])  # "ipopt"
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


@pytest.mark.parametrize("method", ["filtersqp", "funnelsqp", "filterslp"])  # "ipopt"
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
            jac=-2 * np.eye(2),
            lb=np.zeros(2),
            ub=float("Inf"),
        )
    elif c_type == "LinearConstraint":
        constr = LinearConstraint(A=np.eye(2), ub=0.5 * np.zeros(2), lb=-float("Inf"))
    else:
        constr = {"type": "ineq", "fun": lambda x: 1.0 - 2 * x, "jac": -2 * np.eye(2)}
    res = minimize(rosen, jac=rosen_der, x0=np.zeros(2), constraints=(constr,))

    raise ValueError(res)
