# Copyright (c) 2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

import unopy
Inf = float("inf")

# hs015.mod
number_variables = 2
number_constraints = 2

variables_bounds = [(-Inf, 0.5), (-Inf, Inf)]
constraints_bounds = [(1., Inf), (0., Inf)]
	
def evaluate_objective(x):
	return 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2
	
def evaluate_constraints(x, constraints):
	constraints[0] = x[0]*x[1]
	constraints[1] = x[0] + x[1]**2
	
def evaluate_objective_gradient(x, gradient):
	gradient.insert(0, 400.*x[0]*(x[0]**2 - x[0]) + 2.*x[0] - 2.)
	gradient.insert(1, 200.*(x[1] - x[0]**2))

def evaluate_jacobian(x, jacobian):
	# c0
	jacobian.insert(x[1], 0, 0)
	jacobian.insert(x[0], 1, 0)
	# c1
	jacobian.insert(1., 0, 1)
	jacobian.insert(2.*x[1], 1, 1)

def evaluate_hessian(x, objective_multiplier, y, hessian):
	hessian.insert(objective_multiplier*(1200*x[0]**2 - 400.*x[2] + 2.), 0, 0)
	hessian.insert(-400.*objective_multiplier*x[0] + y[0], 0, 1)
	hessian.insert(200.*objective_multiplier + 2.*y[1], 1, 1)

if __name__ == '__main__':
	options = unopy.get_default_options()
	unopy.set_preset(options, "filtersqp")
	constrained_model = (0 < number_constraints)
	uno_solver = unopy.UnoSolver(constrained_model, options)
	result = uno_solver.solve(number_variables, number_constraints,
		evaluate_objective,
		evaluate_constraints,
		evaluate_objective_gradient,
		evaluate_jacobian,
		evaluate_hessian,
		options)
