import pytest
import unopy
from gemseo.algos.design_space import DesignSpace
from gemseo.algos.optimization_problem import OptimizationProblem
from gemseo.algos.opt.factory import OptimizationLibraryFactory
from gemseo.core.mdo_functions.mdo_function import MDOFunction
from gemseo.problems.optimization.power_2 import Power2
from numpy import allclose
from numpy import array

from gemseo_uno.gemseo_uno import UnoOpt
from gemseo_uno.settings.base_uno_settings import UNO_Settings


@pytest.mark.parametrize(
    "method", ["UNO_Filter_SQP", "UNO_Funnel_SQP", "UNO_IPOPT"]
)
def test_power2(method):
    problem = Power2()
    problem.preprocess_functions()
    res = UnoOpt(method).execute(
        problem, settings_model=UNO_Settings(xtol_rel=1e-4, max_iter=50, normalize_design_space=False)
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
        problem, settings_model=UNO_Settings(max_iter=3, normalize_design_space=False)
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
            settings_model=UNO_Settings(max_iter=50, normalize_design_space=False),
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
        settings_model=UNO_Settings(
            max_iter=50,
            xtol_rel=1e-4,
            normalize_design_space=False,
            TR_radius=0.5,
            print_solution=False,
            filter_capacity=10,
        ),
    )
    assert options == {
        "primal_tolerance": 0.0,
        "dual_tolerance": 0.0,
        "TR_radius": 0.5,
        "print_solution": False,
        "filter_capacity": 10,
    }
    assert res.is_feasible


def test_message_and_status():
    """Check that the message is a string and the status an integer."""
    problem = Power2()
    res = UnoOpt("UNO_Filter_SQP").execute(
        problem,
        settings_model=UNO_Settings(max_iter=50, normalize_design_space=False),
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
        settings_model=UNO_Settings(max_iter=100, normalize_design_space=False),
    )
    reference = OptimizationLibraryFactory().execute(
        _create_coupled_problem(), algo_name="SLSQP", max_iter=100, ftol_rel=1e-12
    )
    assert res.is_feasible
    assert allclose(res.x_opt, reference.x_opt, atol=1e-4)
