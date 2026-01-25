import pytest
from gemseo.problems.optimization.power_2 import Power2
from numpy import allclose

from gemseo_uno.gemseo_uno import UnoOpt
from gemseo_uno.settings.base_uno_settings import UNO_Settings


@pytest.mark.parametrize(
    "method", ["UNO_Filter_SQP", "UNO_Funnel_SQP", "UNO_Filter_SLP"]
)
def test_power2(method):
    problem = Power2()
    problem.preprocess_functions()
    res = UnoOpt(method).execute(
        problem, settings_model=UNO_Settings(xtol_rel=1e-4, max_iter=50)
    )
    assert res.is_feasible
    assert allclose(
        res.x_opt, [0.5 ** (1 / 3), 0.5 ** (1 / 3), 0.9 ** (1 / 3)], rtol=1e-3
    )
    assert 3 <= len(problem.database) <= 50
