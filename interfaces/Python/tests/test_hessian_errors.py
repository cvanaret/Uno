# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# evaluation errors raised by the Lagrangian Hessian callbacks

import numpy as np
import pytest
import unopy

Inf = float("inf")


def hs015_positive_convention(hessian, fail_at=None):
	"""hs015 written with the MULTIPLIER_POSITIVE sign convention (L = f + y^T c). The Hessian callback raises at its
	fail_at-th call (1-based). Returns the model and the list of (x, multipliers) received by the Hessian callback."""
	calls = []

	def objective(x):
		return 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2

	def objective_gradient(x, gradient):
		gradient[:] = [400.*x[0]**3 - 400.*x[0]*x[1] + 2.*x[0] - 2., 200.*(x[1] - x[0]**2)]

	def constraints(x, constraint_values):
		constraint_values[:] = [x[0]*x[1], x[0] + x[1]**2]

	def jacobian(x, jacobian_values):
		jacobian_values[:] = [x[1], 1., x[0], 2.*x[1]]

	def hessian_entries(x, objective_multiplier, multipliers):  # lower triangle: (0, 0), (1, 0), (1, 1)
		return (objective_multiplier*(1200.*x[0]**2 - 400.*x[1] + 2.),
				-400.*objective_multiplier*x[0] + multipliers[0],
				200.*objective_multiplier + 2.*multipliers[1])

	def record_call(x, multipliers):
		calls.append((np.array(x), np.array(multipliers)))
		if len(calls) == fail_at:
			raise RuntimeError("injected Hessian evaluation failure")

	def lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values):
		record_call(x, multipliers)
		hessian_values[:] = hessian_entries(x, objective_multiplier, multipliers)

	def lagrangian_hessian_operator(x, evaluate_at_x, objective_multiplier, multipliers, vector, result):
		record_call(x, multipliers)
		h00, h10, h11 = hessian_entries(x, objective_multiplier, multipliers)
		result[:] = [h00*vector[0] + h10*vector[1], h10*vector[0] + h11*vector[1]]

	model = unopy.Model(unopy.PROBLEM_NONLINEAR, 2, unopy.ZERO_BASED_INDEXING)
	model.set_variables_lower_bounds([-Inf, -Inf])
	model.set_variables_upper_bounds([0.5, Inf])
	model.set_objective(unopy.MINIMIZE, objective, objective_gradient)
	model.set_constraints(2, constraints, [1., 0.], [Inf, Inf], 4, [0, 1, 0, 1], [0, 0, 1, 1], jacobian)
	if hessian == "matrix":
		model.set_lagrangian_hessian(3, unopy.LOWER_TRIANGLE, [0, 1, 1], [0, 0, 1], lagrangian_hessian)
	else:
		model.set_lagrangian_hessian_operator(lagrangian_hessian_operator)
	model.set_lagrangian_sign_convention(unopy.MULTIPLIER_POSITIVE)
	model.set_initial_primal_iterate([-2., 1.])
	return model, calls


def filtersqp_solver():
	solver = unopy.UnoSolver()
	solver.set_option("logger", "SILENT")
	solver.set_preset("filtersqp")
	solver.set_option("QP_solver", "BQPD")
	return solver


def test_hessian_operator_failure_preserves_multiplier_signs():
	"""With MULTIPLIER_POSITIVE, the multipliers are flipped for the user callback. A failing Hessian operator must not
	leave Uno's multipliers negated: the calls made at the same point after the failure (the trust-region subproblem is
	solved again at the current iterate) receive the same multipliers."""
	hessian = "operator"
	# reference solve: find the first Hessian call with nonzero multipliers
	model, reference_calls = hs015_positive_convention(hessian)
	filtersqp_solver().optimize(model)
	fail_at = next(index for index, (_, multipliers) in enumerate(reference_calls, start=1) if np.any(multipliers != 0.))

	model, calls = hs015_positive_convention(hessian, fail_at)
	result = filtersqp_solver().optimize(model)
	assert result.optimization_status == unopy.OptimizationStatus.SUCCESS
	assert result.solution_status == unopy.SolutionStatus.FEASIBLE_KKT_POINT
	assert np.allclose(result.primal_solution, [0.5, 2.], atol=1e-4)

	failing_x, failing_multipliers = calls[fail_at - 1]
	assert np.any(failing_multipliers != 0.)
	later_calls_at_same_point = [multipliers for x, multipliers in calls[fail_at:] if np.array_equal(x, failing_x)]
	assert later_calls_at_same_point, "the failing point is never revisited: the test checks nothing"
	for multipliers in later_calls_at_same_point:
		assert np.array_equal(multipliers, failing_multipliers)
