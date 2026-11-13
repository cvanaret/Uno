# GEMSEO interface

[GEMSEO](https://gemseo.readthedocs.io) is a Python library for multidisciplinary design optimization (MDO).
`unopy.gemseo` is a GEMSEO plugin that makes the Uno presets available as GEMSEO optimization algorithms,
so that any GEMSEO optimization problem or MDO scenario can be solved with Uno.

It is part of [unopy](python.md) and its dependencies (GEMSEO 6.3 or later, and SciPy) are installed with the `gemseo` extra:
```bash
pip install "unopy[gemseo]"
```

GEMSEO discovers the plugin automatically once unopy is installed:
the algorithms below are then listed among the optimization algorithms available in GEMSEO.

## Algorithms

| GEMSEO algorithm | Uno preset | Settings class |
|------------------|------------|----------------|
| `UNO_Filter_SQP` | `filtersqp` | `UNO_Filter_SQP_Settings` (`unopy.gemseo.settings.filter_sqp`) |
| `UNO_IPOPT` | `ipopt` | `UNO_IPOPT_Settings` (`unopy.gemseo.settings.ipopt`) |

Both handle equality and inequality constraints and use the gradients provided by GEMSEO.
The functions are only evaluated inside the design space:
Uno can try points slightly outside the bounds (e.g. the IPOPT preset relaxes them),
which are projected onto the bounds before GEMSEO evaluates the functions.
They approximate the Hessian of the Lagrangian with Uno's L-BFGS model.

## Settings

As for any GEMSEO algorithm, the settings of an algorithm are a [Pydantic](https://docs.pydantic.dev) model,
which validates them before the optimization.
The settings of the Uno algorithms are of two kinds:

- the **GEMSEO settings**, common to the gradient-based optimizers of GEMSEO,
  which control the driver: `max_iter`, `max_time`, `ftol_rel`, `ftol_abs`, `xtol_rel`, `xtol_abs`, `stop_crit_n_x`,
  `kkt_tol_abs`, `kkt_tol_rel`, `eq_tolerance`, `ineq_tolerance`, `normalize_design_space`, ...;
- the **Uno options**, one optional field per [Uno option](../options.md) with the same name and type,
  e.g. `TR_radius`, `filter_capacity`, `barrier_initial_parameter`, `linear_solver` or `logger`.
  They default to `None`, which means the default value of Uno (or of the preset);
  only the options set by the user are passed to Uno.
  The options `preset` and `max_iterations` are not available, as the plugin sets them
  (see [Stopping criteria](#stopping-criteria)).

The settings are passed as a settings model, whose class sets the algorithm:

```python
from unopy.gemseo.settings.filter_sqp import UNO_Filter_SQP_Settings

settings = UNO_Filter_SQP_Settings(max_iter=50, xtol_rel=1e-8, TR_radius=0.5, logger="SILENT")
scenario.execute(settings)
```

The settings model gives autocompletion and documentation in an IDE,
and validates the settings:
a misspelled setting or a value of the wrong type raises a `pydantic.ValidationError` before the optimization.

```python
UNO_Filter_SQP_Settings(TR_radiuss=0.5)   # ValidationError: Extra inputs are not permitted
UNO_Filter_SQP_Settings(TR_radius="big")  # ValidationError: Input should be a valid number
```

### Stopping criteria

The stopping criteria are mainly left to GEMSEO (`max_iter`, `max_time`, `ftol_rel`, `xtol_rel`, ...).
To this end, the plugin sets the Uno options `max_iterations` to the largest integer
and `primal_tolerance` and `dual_tolerance` to 1e-15,
so that Uno's own convergence test only stops it once the problem is solved to machine precision.
Setting `primal_tolerance` or `dual_tolerance` overrides these values.
Uno can still stop before GEMSEO with its acceptable test:
after `acceptable_iteration_threshold` (15 by default) consecutive iterates satisfying its looser acceptable tolerances
(`acceptable_tolerance`, 1e-6 by default, and the `acceptable_unscaled_*` options).
This is what stops the IPOPT preset once it no longer makes progress.

## Example: the Power2 problem

The Power2 problem of GEMSEO minimizes $x_1^2 + x_2^2 + x_3^2$
under one equality constraint and two inequality constraints, with $x \in [-1, 1]^3$.
It is an `OptimizationProblem`, solved with `execute_algo`:

```python
from gemseo import execute_algo
from gemseo.problems.optimization.power_2 import Power2

from unopy.gemseo.settings.filter_sqp import UNO_Filter_SQP_Settings

problem = Power2()
settings = UNO_Filter_SQP_Settings(
    max_iter=50,          # GEMSEO setting: maximum number of iterations
    xtol_rel=1e-8,        # GEMSEO setting: relative tolerance on the design variables
    TR_radius=0.5,        # Uno option: initial trust-region radius
    print_solution=False, # Uno option
    logger="SILENT",      # Uno option: no output from Uno
)
result = execute_algo(problem, settings_model=settings)
print(result.x_opt)        # [0.79370053 0.79370053 0.96548938]
print(result.f_opt)        # 2.1920908016810308
print(result.is_feasible)  # True
```

The same problem with the IPOPT preset:

```python
from unopy.gemseo.settings.ipopt import UNO_IPOPT_Settings

problem = Power2()
result = execute_algo(problem, settings_model=UNO_IPOPT_Settings(max_iter=50, logger="SILENT"))
```

## Example: the supersonic business jet (SSBJ)

The SSBJ problem of Sobieski is a classical MDO benchmark of GEMSEO:
four coupled disciplines (propulsion, aerodynamics, structure and mission) and the range of the aircraft to maximize
under structural, aerodynamic and propulsion constraints.
Here it is solved with the MDF formulation, where GEMSEO solves the coupling between the disciplines at each iteration
and computes the coupled derivatives passed to Uno:

```python
from gemseo import create_discipline
from gemseo import create_scenario
from gemseo.problems.mdo.sobieski.core.design_space import SobieskiDesignSpace

from unopy.gemseo.settings.filter_sqp import UNO_Filter_SQP_Settings

disciplines = create_discipline([
    "SobieskiPropulsion",
    "SobieskiAerodynamics",
    "SobieskiMission",
    "SobieskiStructure",
])
scenario = create_scenario(
    disciplines,
    "y_4",  # the range
    SobieskiDesignSpace(),
    maximize_objective=True,
    formulation_name="MDF",
)
for constraint in ["g_1", "g_2", "g_3"]:
    scenario.add_constraint(constraint, constraint_type="ineq")

scenario.execute(UNO_Filter_SQP_Settings(max_iter=20, logger="SILENT"))
result = scenario.optimization_result
print(-result.f_opt)       # 3963.38 (the range is maximized)
print(result.is_feasible)  # True
```

Uno finds the optimal range of 3963.38 in a few iterations.
Another algorithm is used by changing the settings class, e.g. `UNO_IPOPT_Settings`.
