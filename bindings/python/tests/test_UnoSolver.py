# Copyright (c) 2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

import unopy

constrained_model = True
number_variables = 2
number_constraints = 1
	
def evaluate_objective(x):
	return x[0] + x[1]
	
def evaluate_constraints(x, c):
	c[0] = x[0]**2 - x[1]
	
def evaluate_objective_gradient(x, g):
	g.insert(0, 1.)
	g.insert(1, 1.)
	
def evaluate_hessian(x, y, rho, hessian):
	hessian.insert(100., 0, 0)
	hessian.insert(200., 0, 1)
	hessian.insert(300., 1, 1)

if __name__ == '__main__':
	options = unopy.Options.get_default()
	unopy.Options.set_preset(options, "filtersqp")
	solver = unopy.UnoSolver(constrained_model, options)
	solver.solve(number_variables, number_constraints,
		evaluate_objective,
		evaluate_constraints,
		evaluate_objective_gradient,
		evaluate_hessian,
		options)
