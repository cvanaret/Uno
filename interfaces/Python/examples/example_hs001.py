# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# hs001.mod: bound-constrained problem without general constraints. The constraint multipliers passed to the
# Lagrangian Hessian and to the user callbacks are empty arrays (regression test for empty buffers)

import unopy
Inf = float("inf")

def objective(x):
	return 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2

def objective_gradient(x, gradient):
	gradient[:] = [-400.*x[0]*(x[1] - x[0]**2) - 2.*(1. - x[0]),
				   200.*(x[1] - x[0]**2)]

def lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values):
	assert len(multipliers) == 0
	hessian_values[:] = [objective_multiplier*(1200.*x[0]**2 - 400.*x[1] + 2.),
						 -400.*objective_multiplier*x[0],
						 200.*objective_multiplier]

def notify_acceptable_iterate_callback(primals, lower_bound_multipliers, upper_bound_multipliers, constraint_multipliers,
									   objective_multiplier, primal_feasibility_residual, stationarity_residual,
									   complementarity_residual):
	assert len(constraint_multipliers) == 0

def termination_callback(primals, lower_bound_multipliers, upper_bound_multipliers, constraint_multipliers, objective_multiplier,
						 primal_feasibility_residual, stationarity_residual, complementarity_residual):
	assert len(constraint_multipliers) == 0
	return False

if __name__ == '__main__':
	# model creation: 2 variables, x1 >= -1.5, no general constraints
	model = unopy.Model(unopy.PROBLEM_NONLINEAR, 2, unopy.ZERO_BASED_INDEXING)
	model.set_name("hs001")
	model.set_variables_lower_bounds([-Inf, -1.5])
	model.set_variables_upper_bounds([Inf, Inf])
	model.set_objective(unopy.MINIMIZE, objective, objective_gradient)
	model.set_lagrangian_hessian(3, unopy.LOWER_TRIANGLE, [0, 1, 1], [0, 0, 1], lagrangian_hessian)
	model.set_lagrangian_sign_convention(unopy.MULTIPLIER_NEGATIVE)
	model.set_initial_primal_iterate([-2., 1.])

	uno_solver = unopy.UnoSolver()
	uno_solver.set_notify_acceptable_iterate_callback(notify_acceptable_iterate_callback)
	uno_solver.set_termination_callback(termination_callback)
	print("Solving with Uno", unopy.current_uno_version())

	# solution: x* = (1, 1), f(x*) = 0
	for preset, option, solver in [("ipopt", "linear_solver", "MUMPS"), ("filtersqp", "QP_solver", "BQPD")]:
		uno_solver.set_preset(preset)
		uno_solver.set_option(option, solver)
		result = uno_solver.optimize(model)
		print("Method used:", uno_solver.get_method_description())
		print("Optimization status:", result.optimization_status)
		print("Solution status:", result.solution_status)
		print("Objective at solution:", result.solution_objective)
		assert result.optimization_status == unopy.OptimizationStatus.SUCCESS
		assert result.solution_status == unopy.SolutionStatus.FEASIBLE_KKT_POINT
		assert abs(result.solution_objective) <= 1e-6
		assert all(abs(xi - 1.) <= 1e-4 for xi in result.primal_solution)
