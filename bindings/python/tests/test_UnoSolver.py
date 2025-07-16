# Copyright (c) 2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

import unopy
Inf = float("inf")

# hs015.mod
number_variables = 2
number_constraints = 2

variables_lower_bounds = [-Inf, -Inf]
variables_upper_bounds = [0.5, Inf]
constraints_lower_bounds = [1., 0.]
constraints_upper_bounds = [Inf, Inf]
	
def evaluate_objective(x):
	return 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2
	
def evaluate_constraints(x, constraints):
	constraints[0] = x[0]*x[1]
	constraints[1] = x[0] + x[1]**2
	
def evaluate_objective_gradient(x, gradient):
	gradient.insert(0, 400.*x[0]**3 - 400.*x[0]*x[1] + 2.*x[0] - 2.)
	gradient.insert(1, 200.*(x[1] - x[0]**2))

def evaluate_jacobian(x, jacobian):
	# c0
	jacobian[0].insert(0, x[1])
	jacobian[0].insert(1, x[0])
	# c1
	jacobian[1].insert(0, 1.)
	jacobian[1].insert(1, 2.*x[1])

def evaluate_lagrangian_hessian(x, objective_multiplier, y, hessian):
	hessian.insert(0, 0, objective_multiplier*(1200*x[0]**2 - 400.*x[1] + 2.))
	hessian.insert(0, 1, -400.*objective_multiplier*x[0] - y[0])
	hessian.insert(1, 1, 200.*objective_multiplier - 2.*y[1])

primal_initial_point = [-2., 1.]
dual_initial_point = [0., 0.]

if __name__ == '__main__':
	options = unopy.get_default_options()
	unopy.set_preset(options, "filtersqp")
	#options["logger"] = "DEBUG3"
	#options["print_subproblem"] = "yes"
	constrained_model = (0 < number_constraints)
	uno_solver = unopy.UnoSolver(constrained_model, options)
	result = uno_solver.solve(number_variables, number_constraints,
		evaluate_objective,
		evaluate_constraints,
		evaluate_objective_gradient,
		evaluate_jacobian,
		evaluate_lagrangian_hessian,
		variables_lower_bounds,
		variables_upper_bounds,
		constraints_lower_bounds,
		constraints_upper_bounds,
		primal_initial_point,
		dual_initial_point,
		options)
