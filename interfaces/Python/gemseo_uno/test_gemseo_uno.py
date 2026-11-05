import pytest
from gemseo.core.mdo_functions.mdo_function import MDOFunction
from gemseo.problems.optimization.power_2 import Power2
from numpy import allclose

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
