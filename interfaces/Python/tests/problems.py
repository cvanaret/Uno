# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# test problems shared by the unopy tests

from collections import Counter
import unopy

Inf = float("inf")


class CallLog:
	"""Records the callback calls and the sizes of the constraint multipliers they receive. Uno handles exceptions raised
	in the callbacks as evaluation errors, so an assertion inside a callback would not fail a test: record, then assert
	after the solve."""
	def __init__(self):
		self.counts = Counter()
		self.multiplier_sizes = []

	def record(self, name, multipliers=None):
		self.counts[name] += 1
		if multipliers is not None:
			self.multiplier_sizes.append(len(multipliers))


def hs015(hessian="matrix", log=None, triangle=unopy.LOWER_TRIANGLE, convert=list):
	"""hs015: 2 variables, 2 inequality constraints. Solution: x* = (0.5, 2), f(x*) = 306.5.
	hessian: "matrix", "operator" or None. convert: applied to every array argument of the setters."""
	log = log if log is not None else CallLog()

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
				-400.*objective_multiplier*x[0] - multipliers[0],
				200.*objective_multiplier - 2.*multipliers[1])

	def lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values):
		log.record("hessian", multipliers)
		hessian_values[:] = hessian_entries(x, objective_multiplier, multipliers)

	def lagrangian_hessian_operator(x, evaluate_at_x, objective_multiplier, multipliers, vector, result):
		log.record("hessian_operator", multipliers)
		h00, h10, h11 = hessian_entries(x, objective_multiplier, multipliers)
		result[:] = [h00*vector[0] + h10*vector[1], h10*vector[0] + h11*vector[1]]

	model = unopy.Model(unopy.PROBLEM_NONLINEAR, 2, unopy.ZERO_BASED_INDEXING)
	model.set_name("hs015")
	model.set_variables_lower_bounds(convert([-Inf, -Inf]))
	model.set_variables_upper_bounds(convert([0.5, Inf]))
	model.set_objective(unopy.MINIMIZE, objective, objective_gradient)
	model.set_constraints(2, constraints, convert([1., 0.]), convert([Inf, Inf]), 4, convert([0, 1, 0, 1]),
		convert([0, 0, 1, 1]), jacobian)
	if hessian == "matrix":
		rows, columns = [0, 1, 1], [0, 0, 1]
		if triangle == unopy.UPPER_TRIANGLE:
			rows, columns = columns, rows
		model.set_lagrangian_hessian(3, triangle, convert(rows), convert(columns), lagrangian_hessian)
	elif hessian == "operator":
		model.set_lagrangian_hessian_operator(lagrangian_hessian_operator)
	model.set_lagrangian_sign_convention(unopy.MULTIPLIER_NEGATIVE)
	model.set_initial_primal_iterate(convert([-2., 1.]))
	return model


def hs001(log=None):
	"""hs001: 2 variables, bound constraint x[1] >= -1.5, no general constraints. Solution: x* = (1, 1), f(x*) = 0"""
	log = log if log is not None else CallLog()

	def objective(x):
		return 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2

	def objective_gradient(x, gradient):
		gradient[:] = [-400.*x[0]*(x[1] - x[0]**2) - 2.*(1. - x[0]), 200.*(x[1] - x[0]**2)]

	def lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values):
		log.record("hessian", multipliers)
		hessian_values[:] = [objective_multiplier*(1200.*x[0]**2 - 400.*x[1] + 2.), -400.*objective_multiplier*x[0],
							 200.*objective_multiplier]

	model = unopy.Model(unopy.PROBLEM_NONLINEAR, 2, unopy.ZERO_BASED_INDEXING)
	model.set_name("hs001")
	model.set_variables_lower_bounds([-Inf, -1.5])
	model.set_variables_upper_bounds([Inf, Inf])
	model.set_objective(unopy.MINIMIZE, objective, objective_gradient)
	model.set_lagrangian_hessian(3, unopy.LOWER_TRIANGLE, [0, 1, 1], [0, 0, 1], lagrangian_hessian)
	model.set_initial_primal_iterate([-2., 1.])
	return model


def record_solver_callbacks(solver, log):
	"""Records the calls of the notify and termination callbacks (constraint multipliers: 4th argument)"""
	solver.set_notify_acceptable_iterate_callback(lambda *arguments: log.record("notify", arguments[3]))
	solver.set_termination_callback(lambda *arguments: log.record("termination", arguments[3]) or False)
