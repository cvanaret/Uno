# Copyright (c) 2021-2026 IRT Saint Exupéry, https://www.irt-saintexupery.com
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""Tests of the GEMSEO plugin unopy.gemseo."""

from __future__ import annotations

import re
from types import NoneType
from typing import get_args

import pytest

pytest.importorskip("gemseo")

import unopy  # noqa: E402
import unopy.gemseo.uno as uno_module  # noqa: E402
from gemseo.algos.design_space import DesignSpace  # noqa: E402
from gemseo.algos.opt.base_gradient_based_algorithm_settings import (  # noqa: E402
    BaseGradientBasedAlgorithmSettings,
)
from gemseo.algos.opt.base_optimizer_settings import BaseOptimizerSettings  # noqa: E402
from gemseo.algos.opt.factory import OptimizationLibraryFactory  # noqa: E402
from gemseo.algos.optimization_problem import OptimizationProblem  # noqa: E402
from gemseo.core.mdo_functions.mdo_function import MDOFunction  # noqa: E402
from gemseo.problems.optimization.power_2 import Power2  # noqa: E402
from gemseo.utils.constants import C_LONG_MAX  # noqa: E402
from numpy import allclose  # noqa: E402
from numpy import array  # noqa: E402
from numpy import identity  # noqa: E402
from numpy import vstack  # noqa: E402
from numpy.testing import assert_array_equal  # noqa: E402
from scipy.sparse import csr_array  # noqa: E402
from scipy.sparse import issparse  # noqa: E402

from unopy.gemseo.settings.base_uno_settings import BaseUnoSettings  # noqa: E402
from unopy.gemseo.settings.filter_slp import UNO_Filter_SLP_Settings  # noqa: E402
from unopy.gemseo.settings.filter_sqp import UNO_Filter_SQP_Settings  # noqa: E402
from unopy.gemseo.settings.funnel_sqp import UNO_Funnel_SQP_Settings  # noqa: E402
from unopy.gemseo.settings.ipopt import UNO_IPOPT_Settings  # noqa: E402
from unopy.gemseo.uno import UnoOpt  # noqa: E402
from unopy.gemseo.uno import _concatenate_constraints  # noqa: E402
from unopy.scipy import minimize  # noqa: E402

SETTINGS = {
    "UNO_Filter_SQP": UNO_Filter_SQP_Settings,
    "UNO_Filter_SLP": UNO_Filter_SLP_Settings,
    "UNO_Funnel_SQP": UNO_Funnel_SQP_Settings,
    "UNO_IPOPT": UNO_IPOPT_Settings,
}
"""The settings classes of the Uno algorithms."""


@pytest.mark.parametrize(
    "method", ["UNO_Filter_SQP", "UNO_Funnel_SQP", "UNO_IPOPT"]
)
def test_power2(method):
    problem = Power2()
    problem.preprocess_functions()
    res = UnoOpt(method).execute(
        problem, settings_model=SETTINGS[method](xtol_rel=1e-4, max_iter=50, normalize_design_space=False)
    )
    assert res.is_feasible
    assert allclose(
        res.x_opt, [0.5 ** (1 / 3), 0.5 ** (1 / 3), 0.9 ** (1 / 3)], rtol=1e-3
    )
    assert 3 <= len(problem.database) <= 50


@pytest.mark.parametrize(
    "method", ["UNO_Filter_SQP", "UNO_Funnel_SQP", "UNO_IPOPT"]
)
def test_power2_max_iter(method):
    """Check that GEMSEO stops Uno when the maximum number of iterations is reached."""
    problem = Power2()
    res = UnoOpt(method).execute(
        problem, settings_model=SETTINGS[method](max_iter=3, normalize_design_space=False)
    )
    assert len(problem.database) == 3
    assert res.message == (
        "Maximum number of iterations reached. GEMSEO stopped the driver."
    )


class UserError(Exception):
    """An error raised by a user function."""


@pytest.mark.parametrize(
    "method", ["UNO_Filter_SQP", "UNO_Funnel_SQP", "UNO_IPOPT"]
)
def test_power2_user_exception(method):
    """Check that an exception raised by a function is propagated."""
    problem = Power2()
    calls = []

    def failing_function(x):
        calls.append(x)
        if len(calls) > 2:
            raise UserError("Failure.")
        return (x**2).sum()

    problem.objective = MDOFunction(
        failing_function, "pow2", jac=lambda x: 2 * x
    )
    with pytest.raises(UserError, match="Failure."):
        UnoOpt(method).execute(
            problem,
            settings_model=SETTINGS[method](max_iter=50, normalize_design_space=False),
        )
    assert len(calls) == 3


def test_settings_forwarded_to_uno(monkeypatch):
    """Check that only the Uno settings set by the user are passed to Uno."""
    options = {}
    set_option = unopy.UnoSolver.set_option

    def spy(solver, name, value):
        options[name] = value
        set_option(solver, name, value)

    monkeypatch.setattr(unopy.UnoSolver, "set_option", spy)
    problem = Power2()
    res = UnoOpt("UNO_Filter_SQP").execute(
        problem,
        settings_model=UNO_Filter_SQP_Settings(
            max_iter=50,
            xtol_rel=1e-4,
            normalize_design_space=False,
            TR_radius=0.5,
            print_solution=False,
            filter_capacity=10,
            SOC_max_iterations=3,
            barrier_small_infeasibility_factor=0.5,
        ),
    )
    assert options == {
        "primal_tolerance": 0.0,
        "dual_tolerance": 0.0,
        "TR_radius": 0.5,
        "print_solution": False,
        "filter_capacity": 10,
        "SOC_max_iterations": 3,
        "barrier_small_infeasibility_factor": 0.5,
        "max_iterations": C_LONG_MAX,
    }
    assert res.is_feasible


def test_uno_settings_are_uno_options():
    """Check that the Uno settings are Uno options with the right types.

    The Uno options and their types are read from the documentation of minimize,
    which lists the options of uno/options/Options.cpp.
    """
    documented_options = {
        name: {"float": float, "int": int, "bool": bool, "str": str}[type_name]
        for name, type_name in re.findall(
            r'^\s*"(\w+)"[^:\n]*: (\w+)\s*$', minimize.__doc__, flags=re.MULTILINE
        )
    }
    base_fields = set(BaseOptimizerSettings.model_fields) | set(
        BaseGradientBasedAlgorithmSettings.model_fields
    )
    uno_fields = {
        name: field
        for name, field in BaseUnoSettings.model_fields.items()
        if name not in base_fields
    }
    assert uno_fields
    for name, field in uno_fields.items():
        assert name in documented_options, name
        assert set(get_args(field.annotation)) == {documented_options[name], NoneType}
        assert field.default is None
    # All the options except those set by the wrapper are exposed.
    assert set(documented_options) - set(uno_fields) == {"max_iterations"}


def test_max_iter_above_uno_default(monkeypatch):
    """Check that Uno's iteration limit does not stop it before GEMSEO's max_iter."""
    options = {}
    set_option = unopy.UnoSolver.set_option

    def spy(solver, name, value):
        options[name] = value
        set_option(solver, name, value)

    monkeypatch.setattr(unopy.UnoSolver, "set_option", spy)
    # Uno's default maximum number of iterations is 2000.
    UnoOpt("UNO_Filter_SQP").execute(
        Power2(),
        settings_model=UNO_Filter_SQP_Settings(max_iter=5000, normalize_design_space=False),
    )
    assert options["max_iterations"] >= 5000


def test_message_and_status():
    """Check that the message is a string and the status an integer."""
    problem = Power2()
    res = UnoOpt("UNO_Filter_SQP").execute(
        problem,
        settings_model=UNO_Filter_SQP_Settings(max_iter=50, normalize_design_space=False),
    )
    assert isinstance(res.message, str)
    assert type(res.status) is int


def _create_coupled_problem():
    """Create a problem with several constraints with non-diagonal Jacobians."""
    design_space = DesignSpace()
    design_space.add_variable("x", 3, lower_bound=-10.0, upper_bound=10.0, value=0.0)
    problem = OptimizationProblem(design_space)
    problem.objective = MDOFunction(
        lambda x: ((x - 2.0) ** 2).sum(), "f", jac=lambda x: 2 * (x - 2.0)
    )
    a_1 = array([1.0, 2.0, 3.0])
    a_2 = array([1.0, 0.0, -4.0])
    a_3 = array([0.5, 1.0, 1.0])
    for name, a, value in (("c_1", a_1, 2.0), ("c_2", a_2, 1.0), ("c_3", a_3, 1.0)):
        problem.add_constraint(
            MDOFunction(
                lambda x, a=a, value=value: array([a @ x - value]),
                name,
                jac=lambda x, a=a: a.reshape(1, -1),
            ),
            constraint_type=MDOFunction.ConstraintType.INEQ,
        )
    return problem


@pytest.mark.parametrize(
    "method", ["UNO_Filter_SQP", "UNO_Funnel_SQP", "UNO_IPOPT"]
)
def test_several_constraints(method):
    """Check a problem with several constraints having non-diagonal Jacobians."""
    problem = _create_coupled_problem()
    res = UnoOpt(method).execute(
        problem,
        settings_model=SETTINGS[method](max_iter=100, normalize_design_space=False),
    )
    reference = OptimizationLibraryFactory().execute(
        _create_coupled_problem(), algo_name="SLSQP", max_iter=100, ftol_rel=1e-12
    )
    assert res.is_feasible
    assert allclose(res.x_opt, reference.x_opt, atol=1e-4)


@pytest.mark.parametrize("algo_name", sorted(SETTINGS))
def test_settings_classes(algo_name):
    """Check that each algorithm has its own settings class targeting it."""
    settings_class = UnoOpt.ALGORITHM_INFOS[algo_name].Settings
    assert settings_class is SETTINGS[algo_name]
    assert issubclass(settings_class, BaseUnoSettings)
    assert settings_class._TARGET_CLASS_NAME == algo_name
    assert set(settings_class.model_fields) == set(BaseUnoSettings.model_fields)


@pytest.mark.parametrize("algo_name", sorted(SETTINGS))
def test_plugin_discovery(algo_name):
    """Check that GEMSEO finds the Uno algorithms through the unopy entry point."""
    assert OptimizationLibraryFactory().is_available(algo_name)


def _create_mixed_problem(calls=None):
    """Create a problem with several equality and inequality constraints.

    Some constraints are vector-valued and one is a positive inequality constraint.

    Args:
        calls: The names of the constraints at each call, if any.
    """
    design_space = DesignSpace()
    design_space.add_variable("x", 5, lower_bound=-10.0, upper_bound=10.0, value=0.0)
    problem = OptimizationProblem(design_space)
    problem.objective = MDOFunction(
        lambda x: ((x - 2.0) ** 2).sum(), "f", jac=lambda x: 2 * (x - 2.0)
    )

    def counted(name, function):
        def counted_function(x):
            if calls is not None:
                calls.append(name)
            return function(x)

        return counted_function

    eye = identity(5)
    constraints = (
        # Equality: x0 + x1 = 1 and x2 = x3.
        (
            "eq_vector",
            lambda x: array([x[0] + x[1] - 1.0, x[2] - x[3]]),
            lambda x: vstack([eye[0] + eye[1], eye[2] - eye[3]]),
            MDOFunction.ConstraintType.EQ,
            False,
        ),
        # Inequality: x0 <= 0.3 and x2 + x3 <= 1.
        (
            "ineq_vector",
            lambda x: array([x[0] - 0.3, x[2] + x[3] - 1.0]),
            lambda x: vstack([eye[0], eye[2] + eye[3]]),
            MDOFunction.ConstraintType.INEQ,
            False,
        ),
        # Equality: x4 = 0.5 x0 + 1.
        (
            "eq_scalar",
            lambda x: array([x[4] - 0.5 * x[0] - 1.0]),
            lambda x: (eye[4] - 0.5 * eye[0]).reshape(1, -1),
            MDOFunction.ConstraintType.EQ,
            False,
        ),
        # Positive inequality: 0.5 - x1^2 >= 0.
        (
            "ineq_positive",
            lambda x: array([0.5 - x[1] ** 2]),
            lambda x: (-2 * x[1] * eye[1]).reshape(1, -1),
            MDOFunction.ConstraintType.INEQ,
            True,
        ),
    )
    for name, function, jacobian, constraint_type, positive in constraints:
        problem.add_constraint(
            MDOFunction(counted(name, function), name, jac=jacobian),
            constraint_type=constraint_type,
            positive=positive,
        )
    return problem


@pytest.mark.parametrize(
    "method", ["UNO_Filter_SQP", "UNO_Funnel_SQP", "UNO_IPOPT"]
)
def test_merged_constraints(monkeypatch, method):
    """Check that the constraints are merged into one block per type."""
    scipy_constraints = []

    def spy(*args, **kwargs):
        scipy_constraints.append(kwargs["constraints"])
        return minimize(*args, **kwargs)

    monkeypatch.setattr(uno_module, "minimize", spy)
    calls = []
    problem = _create_mixed_problem(calls)
    res = UnoOpt(method).execute(
        problem,
        settings_model=SETTINGS[method](
            max_iter=100, xtol_rel=1e-10, normalize_design_space=False
        ),
    )
    reference = OptimizationLibraryFactory().execute(
        _create_mixed_problem(), algo_name="SLSQP", max_iter=100, ftol_rel=1e-12
    )
    assert res.is_feasible
    assert allclose(res.x_opt, reference.x_opt, atol=1e-4)
    # One equality and one inequality constraint for the scipy layer.
    assert len(scipy_constraints) == 1
    assert [constraint["type"] for constraint in scipy_constraints[0]] == [
        "eq",
        "ineq",
    ]
    # Each constraint function is evaluated at most once per point of the database.
    for name in ("eq_vector", "ineq_vector", "eq_scalar", "ineq_positive"):
        assert 0 < calls.count(name) <= len(problem.database)


def test_concatenate_constraints():
    """Check the concatenation of the values and the Jacobians of constraints."""
    x = array([1.0, 2.0, 3.0])
    dense = MDOFunction(
        lambda x: x[0] * x[1], "dense", jac=lambda x: array([x[1], x[0], 0.0])
    )
    sparse = MDOFunction(
        lambda x: array([x[2], 2 * x[0]]),
        "sparse",
        jac=lambda x: csr_array(array([[0.0, 0.0, 1.0], [2.0, 0.0, 0.0]])),
    )
    expected_jacobian = array([[2.0, 1.0, 0.0], [0.0, 0.0, 1.0], [2.0, 0.0, 0.0]])

    constraint = _concatenate_constraints([dense, sparse], "ineq")
    assert constraint["type"] == "ineq"
    assert_array_equal(constraint["fun"](x), [2.0, 3.0, 2.0])
    jacobian = constraint["jac"](x)
    assert issparse(jacobian)
    assert_array_equal(jacobian.toarray(), expected_jacobian)

    # Dense Jacobians only: a dense Jacobian.
    constraint = _concatenate_constraints([dense, dense], "eq")
    jacobian = constraint["jac"](x)
    assert not issparse(jacobian)
    assert_array_equal(jacobian, expected_jacobian[[0, 0]])
