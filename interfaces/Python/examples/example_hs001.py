# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# hs001.mod: bound-constrained problem without general constraints. Solution: x* = (1, 1), f(x*) = 0

import unopy
Inf = float("inf")

def objective(x):
	return 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2

def objective_gradient(x, gradient):
	gradient[:] = [-400.*x[0]*(x[1] - x[0]**2) - 2.*(1. - x[0]),
				   200.*(x[1] - x[0]**2)]

def lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values):
	# lower triangle of the Hessian of the objective (no general constraints: the multipliers array is empty)
	hessian_values[:] = [objective_multiplier*(1200.*x[0]**2 - 400.*x[1] + 2.),
						 -400.*objective_multiplier*x[0],
						 200.*objective_multiplier]

if __name__ == '__main__':
	# model creation: 2 variables, x[1] >= -1.5, no general constraints
	model = unopy.Model(unopy.PROBLEM_NONLINEAR, 2, unopy.ZERO_BASED_INDEXING)
	model.set_name("hs001")
	model.set_variables_lower_bounds([-Inf, -1.5])
	model.set_variables_upper_bounds([Inf, Inf])
	model.set_objective(unopy.MINIMIZE, objective, objective_gradient)
	model.set_lagrangian_hessian(3, unopy.LOWER_TRIANGLE, [0, 1, 1], [0, 0, 1], lagrangian_hessian)
	model.set_lagrangian_sign_convention(unopy.MULTIPLIER_NEGATIVE)
	model.set_initial_primal_iterate([-2., 1.])

	uno_solver = unopy.UnoSolver()
	print("Solving with Uno", unopy.current_uno_version())

	# solve with an interior-point method and with an SQP method
	for preset, option, solver in [("ipopt", "linear_solver", "MUMPS"), ("filtersqp", "QP_solver", "BQPD")]:
		uno_solver.set_preset(preset)
		uno_solver.set_option(option, solver)
		result = uno_solver.optimize(model)
		print("Method used:", uno_solver.get_method_description())
		print("Optimization status:", result.optimization_status)
		print("Solution status:", result.solution_status)
		print("Objective at solution:", result.solution_objective)
		print("Primal solution:", result.primal_solution)
