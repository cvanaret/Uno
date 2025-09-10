# Copyright (c) 2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

import unopy
Inf = float("inf")

# hs015.mod

def objective(x, objective_value):
	objective_value[0] = 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2
	return 0
	
def constraints(x, constraint_values):
	constraint_values[0] = x[0]*x[1]
	constraint_values[1] = x[0] + x[1]**2
	return 0
	
def objective_gradient(x, gradient):
	gradient[0] = 400.*x[0]**3 - 400.*x[0]*x[1] + 2.*x[0] - 2.
	gradient[1] = 200.*(x[1] - x[0]**2)
	return 0

def constraint_jacobian(x, jacobian_values):
	jacobian_values[0] = x[1]
	jacobian_values[1] = 1.
	jacobian_values[2] = x[0]
	jacobian_values[3] = 2.*x[1]
	return 0

def lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values):
	hessian_values[0] = objective_multiplier*(1200*x[0]**2 - 400.*x[1] + 2.)
	hessian_values[1] = -400.*objective_multiplier*x[0] - multipliers[0]
	hessian_values[2] = 200.*objective_multiplier - 2.*multipliers[1]
	return 0

if __name__ == '__main__':
	# model creation
	base_indexing = unopy.ZERO_BASED_INDEXING
	# variables
	number_variables = 2
	variables_lower_bounds = [-Inf, -Inf]
	variables_upper_bounds = [0.5, Inf]
	# objective
	optimization_sense = unopy.MINIMIZE
	# constraints
	number_constraints = 2
	constraints_lower_bounds = [1., 0.]
	constraints_upper_bounds = [Inf, Inf]
	number_jacobian_nonzeros = 4
	jacobian_row_indices = [0, 1, 0, 1]
	jacobian_column_indices = [0, 0, 1, 1]
	# Hessian
	number_hessian_nonzeros = 3
	hessian_triangular_part = unopy.LOWER_TRIANGLE
	hessian_row_indices = [0, 1, 1]
	hessian_column_indices = [0, 0, 1]
	lagrangian_sign_convention = unopy.MULTIPLIER_NEGATIVE
	# initial point
	x0 = [-2., 1.]

	model = unopy.Model(unopy.PROBLEM_NONLINEAR, number_variables, variables_lower_bounds, variables_upper_bounds,
		base_indexing)
	model.set_objective(optimization_sense, objective, objective_gradient)
	model.set_constraints(number_constraints, constraints, constraints_lower_bounds, constraints_upper_bounds,
	  number_jacobian_nonzeros, jacobian_row_indices, jacobian_column_indices, constraint_jacobian)
	model.set_lagrangian_hessian(number_hessian_nonzeros, hessian_triangular_part, hessian_row_indices,
		hessian_column_indices, lagrangian_hessian, lagrangian_sign_convention)
	model.set_initial_primal_iterate(x0)
	
	uno_solver = unopy.UnoSolver()
	uno_solver.set_preset("filtersqp")
	#uno_solver.set_option("print_subproblem", "yes")
	#uno_solver.set_option("print_solution", "yes")
	#uno_solver.set_option("logger", "DEBUG3")
	
	result = uno_solver.optimize(model)
	