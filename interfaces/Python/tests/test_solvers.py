# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# end-to-end solves with the presets and subproblem solvers shipped in the wheels

import numpy as np
import pytest
import unopy
from problems import CallLog, hs001, hs015, record_solver_callbacks

SUCCESS = unopy.OptimizationStatus.SUCCESS
FEASIBLE_KKT_POINT = unopy.SolutionStatus.FEASIBLE_KKT_POINT

# (preset, subproblem solver option, subproblem solver, Hessian representation in the model,
#  expected method description, Hessian callback expected to be called)
HS015_METHODS = [
	pytest.param("filtersqp", "QP_solver", "BQPD", None, "L-BFGS Hessian", None, id="filtersqp-BQPD-LBFGS"),
	pytest.param("filtersqp", "QP_solver", "BQPD", "matrix", "exact Hessian", "hessian", id="filtersqp-BQPD-matrix"),
	pytest.param("filtersqp", "QP_solver", "BQPD", "operator", "exact Hessian", "hessian_operator",
		id="filtersqp-BQPD-operator"),
	pytest.param("ipopt", "linear_solver", "MUMPS", "matrix", "interior-point", "hessian", id="ipopt-MUMPS"),
	pytest.param("filterslp", "LP_solver", "HiGHS", "matrix", "zero Hessian", None, id="filterslp-HiGHS"),
]


@pytest.mark.parametrize("preset, option, subproblem_solver, hessian, method, hessian_callback", HS015_METHODS)
def test_hs015(solver, preset, option, subproblem_solver, hessian, method, hessian_callback):
	log = CallLog()
	solver.set_preset(preset)
	solver.set_option(option, subproblem_solver)
	result = solver.optimize(hs015(hessian, log))
	assert method in solver.get_method_description()
	assert result.optimization_status == SUCCESS
	assert result.solution_status == FEASIBLE_KKT_POINT
	assert result.solution_objective == pytest.approx(306.5, abs=1e-4)
	assert np.allclose(result.primal_solution, [0.5, 2.], atol=1e-4)
	# the Hessian is evaluated through the expected representation only
	for callback in ("hessian", "hessian_operator"):
		assert (log.counts[callback] > 0) == (callback == hessian_callback), callback


def test_solver_reuse(solver):
	"""A solver reused across models and presets gives the same results as a fresh solver"""
	solver.set_preset("filtersqp")
	solver.set_option("QP_solver", "BQPD")
	solver.optimize(hs015(hessian=None))
	solver.set_preset("ipopt")
	solver.set_option("linear_solver", "MUMPS")
	reused = solver.optimize(hs015())

	fresh_solver = unopy.UnoSolver()
	fresh_solver.set_option("logger", "SILENT")
	fresh_solver.set_preset("ipopt")
	fresh_solver.set_option("linear_solver", "MUMPS")
	fresh = fresh_solver.optimize(hs015())
	assert reused.number_iterations == fresh.number_iterations
	assert np.array_equal(reused.primal_solution, fresh.primal_solution)


@pytest.mark.parametrize("preset, option, subproblem_solver", [
	("ipopt", "linear_solver", "MUMPS"), ("filtersqp", "QP_solver", "BQPD")])
def test_hs001_without_constraints(solver, preset, option, subproblem_solver):
	"""Without general constraints, the Hessian and the solver callbacks receive empty constraint multipliers"""
	log = CallLog()
	record_solver_callbacks(solver, log)
	solver.set_preset(preset)
	solver.set_option(option, subproblem_solver)
	result = solver.optimize(hs001(log))
	assert result.optimization_status == SUCCESS
	assert result.solution_status == FEASIBLE_KKT_POINT
	assert abs(result.solution_objective) <= 1e-6
	assert np.allclose(result.primal_solution, [1., 1.], atol=1e-4)
	assert len(result.constraint_dual_solution) == 0
	assert log.counts["hessian"] > 0 and log.counts["notify"] > 0
	assert log.multiplier_sizes and all(size == 0 for size in log.multiplier_sizes)
