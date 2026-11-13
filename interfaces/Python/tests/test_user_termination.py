# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# a callback raises unopy.UserTermination to stop Uno immediately

import numpy as np
import pytest
import unopy
from problems import hs015

USER_TERMINATION = unopy.OptimizationStatus.USER_TERMINATION

PRESETS = [
	pytest.param("filtersqp", "QP_solver", "BQPD", id="filtersqp-BQPD"),
	pytest.param("ipopt", "linear_solver", "MUMPS", id="ipopt-MUMPS"),
]


class StopRequest(unopy.UserTermination):
	"""A subclass carrying the user's own data (e.g. the exception to re-raise after optimize)"""
	def __init__(self, reason):
		super().__init__(reason)
		self.reason = reason


def make_solver(preset, option, subproblem_solver):
	solver = unopy.UnoSolver()
	solver.set_option("logger", "SILENT")
	solver.set_preset(preset)
	solver.set_option(option, subproblem_solver)
	return solver


def interruptible_hs015(fail_at, exception):
	"""hs015 whose objective raises exception from its fail_at-th call on (1-based, None: never).
	Returns the model and the number of objective calls (a list with one element)."""
	model = hs015()
	calls = [0]
	objective = None

	def capture_objective(x):
		return 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2

	def objective_gradient(x, gradient):
		gradient[:] = [400.*x[0]**3 - 400.*x[0]*x[1] + 2.*x[0] - 2., 200.*(x[1] - x[0]**2)]

	def interruptible_objective(x):
		calls[0] += 1
		if fail_at is not None and fail_at <= calls[0]:
			raise exception
		return capture_objective(x)

	model.set_objective(unopy.MINIMIZE, interruptible_objective, objective_gradient)
	return model, calls


def record_accepted_iterates(solver, calls):
	"""Records the accepted primals (of the reformulated problem: the original variables come first) and the number of
	objective calls when they are accepted"""
	accepted = []
	solver.set_notify_acceptable_iterate_callback(
		lambda primals, *arguments: accepted.append((np.array(primals), calls[0])))
	return accepted


@pytest.mark.parametrize("preset, option, subproblem_solver", PRESETS)
@pytest.mark.parametrize("exception", [unopy.UserTermination("stop"), StopRequest("stop")], ids=["base", "subclass"])
def test_user_termination_in_model_callback(preset, option, subproblem_solver, exception):
	"""Raising unopy.UserTermination in a model callback stops Uno at once, at the last accepted iterate"""
	# reference run
	model, calls = interruptible_hs015(None, None)
	solver = make_solver(preset, option, subproblem_solver)
	reference_accepted = record_accepted_iterates(solver, calls)
	reference = solver.optimize(model)
	assert reference.optimization_status == unopy.OptimizationStatus.SUCCESS
	assert len(reference_accepted) >= 3

	# the first objective evaluation after the second accepted iterate requests the termination
	fail_at = reference_accepted[1][1] + 1
	model, calls = interruptible_hs015(fail_at, exception)
	solver = make_solver(preset, option, subproblem_solver)
	accepted = record_accepted_iterates(solver, calls)
	result = solver.optimize(model)
	assert result.optimization_status == USER_TERMINATION
	# no evaluation after the termination request (no backtracking, no dummy evaluations)
	assert calls[0] == fail_at
	assert len(accepted) == 2
	assert np.array_equal(result.primal_solution, accepted[-1][0][:2])
	assert np.array_equal(accepted[-1][0], reference_accepted[1][0])
	assert np.isfinite(result.solution_objective)


@pytest.mark.parametrize("preset, option, subproblem_solver", PRESETS)
def test_other_exception_is_an_evaluation_error(preset, option, subproblem_solver):
	"""Contrast: any other exception raised in a model callback is an evaluation error, handled by backtracking"""
	model, calls = interruptible_hs015(None, None)
	solver = make_solver(preset, option, subproblem_solver)
	reference_accepted = record_accepted_iterates(solver, calls)
	solver.optimize(model)

	fail_at = reference_accepted[1][1] + 1
	failures = []

	def failing_once_objective(x):
		calls[0] += 1
		if calls[0] == fail_at:
			failures.append(calls[0])
			raise RuntimeError("evaluation failure")
		return 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2

	model, calls = interruptible_hs015(None, None)
	model.set_objective(unopy.MINIMIZE, failing_once_objective,
		lambda x, g: g.__setitem__(slice(None), [400.*x[0]**3 - 400.*x[0]*x[1] + 2.*x[0] - 2., 200.*(x[1] - x[0]**2)]))
	result = make_solver(preset, option, subproblem_solver).optimize(model)
	assert failures == [fail_at]
	assert result.optimization_status == unopy.OptimizationStatus.SUCCESS
	assert calls[0] > fail_at


@pytest.mark.parametrize("preset, option, subproblem_solver", PRESETS)
def test_user_termination_in_solver_callback(preset, option, subproblem_solver):
	"""Raising unopy.UserTermination in the notify callback also stops Uno at once"""
	notifications = []

	def notify(primals, *arguments):
		notifications.append(np.array(primals))
		if len(notifications) == 2:
			raise unopy.UserTermination()

	solver = make_solver(preset, option, subproblem_solver)
	solver.set_notify_acceptable_iterate_callback(notify)
	result = solver.optimize(hs015())
	assert result.optimization_status == USER_TERMINATION
	assert len(notifications) == 2
	# the iterate being notified was not accepted yet
	assert np.array_equal(result.primal_solution, notifications[0][:2])


def test_user_termination_at_initial_point():
	model, calls = interruptible_hs015(1, unopy.UserTermination())
	result = make_solver("ipopt", "linear_solver", "MUMPS").optimize(model)
	assert result.optimization_status == USER_TERMINATION
	assert result.number_iterations == 0
	# the objective is not called again (e.g. to evaluate the final iterate)
	assert calls[0] == 1
