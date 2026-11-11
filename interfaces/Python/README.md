## unopy

`unopy` is the Python interface for [Uno](https://github.com/cvanaret/Uno), a modern and modular solver for nonlinearly constrained optimization.
Uno unifies Lagrange-Newton (**SQP** and **interior-point**) methods by breaking them down into a set of common building blocks such as constraint reformulation, step computation, and globalization.

## Installation

`unopy` is a registered Python package that can be installed via pip:

```bash
pip install unopy
```

unopy allows you to solve an optimization model described by callback functions.
An example of the Hock-Schittkowski model [hs015](https://vanderbei.princeton.edu/ampl/nlmodels/hs/hs015.mod) is available in the file [example_hs015.py](https://github.com/cvanaret/Uno/blob/main/interfaces/Python/examples/example_hs015.py).

## SciPy and GEMSEO interfaces

unopy ships two optional interfaces, whose dependencies are installed with extras:

```bash
pip install "unopy[scipy]"   # scipy.optimize.minimize-like interface
pip install "unopy[gemseo]"  # GEMSEO plugin
```

`unopy.scipy.minimize` follows the API of `scipy.optimize.minimize`, with the Uno presets as methods (`"filtersqp"` and `"ipopt"`):

```python
from scipy.optimize import rosen, rosen_der
from unopy.scipy import minimize

result = minimize(rosen, [1.3, 0.7, 0.8], jac=rosen_der, method="filtersqp")
```

`unopy.gemseo` is a [GEMSEO](https://gemseo.readthedocs.io) plugin discovered automatically by GEMSEO: the algorithms `UNO_Filter_SQP`, `UNO_Filter_SLP`, `UNO_Funnel_SQP` and `UNO_IPOPT` can be used as any GEMSEO optimization algorithm, with their settings models `UNO_Filter_SQP_Settings`, `UNO_Filter_SLP_Settings`, `UNO_Funnel_SQP_Settings` and `UNO_IPOPT_Settings`:

```python
from unopy.gemseo.settings.filter_sqp import UNO_Filter_SQP_Settings

scenario.execute(UNO_Filter_SQP_Settings(max_iter=100))
```

## Getting started

To get started with Uno, check out the [official documentation](https://unosolver.readthedocs.io/en/latest/interfaces/python/).

## Running the tests

The tests are in [interfaces/Python/tests](https://github.com/cvanaret/Uno/tree/main/interfaces/Python/tests) and use pytest.
From a clone of the repository, download the dependencies, install unopy with the `test` extra and run pytest from the root directory:

```bash
bash dependencies/scripts/download_dependencies.sh
pip install ".[test]"
pytest
```

The tests require the solvers shipped in the wheels: BQPD, HiGHS and MUMPS.
The tests of `unopy.scipy` and `unopy.gemseo` are skipped when SciPy or GEMSEO is not installed; to run them, install the extras, e.g. `pip install ".[test,gemseo]"`.
