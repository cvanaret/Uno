# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# The tests require an unopy build with the solvers shipped in the wheels: BQPD, HiGHS and MUMPS.
# No __init__.py in this directory: otherwise pytest would put interfaces/Python on sys.path, and the source
# directory interfaces/Python/unopy (without the compiled module) would shadow the installed unopy.

import pytest
import unopy


@pytest.fixture
def solver():
	uno_solver = unopy.UnoSolver()
	uno_solver.set_option("logger", "SILENT")
	return uno_solver
