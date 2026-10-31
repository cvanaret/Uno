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
