# Copyright (c) 2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

import unopy
Inf = float("inf")

# hs015.mod

def objective(number_variables, x, objective_value, user_data):
	# objective definition depends on a parameter "offset" passed via user_data
	objective_value[0] = 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2 + user_data["offset"]
	return 0
	
def constraints(number_variables, number_constraints, x, constraint_values, user_data):
	constraint_values[0] = x[0]*x[1]
	constraint_values[1] = x[0] + x[1]**2
	return 0
	
def objective_gradient(number_variables, x, gradient, user_data):
	gradient[0] = 400.*x[0]**3 - 400.*x[0]*x[1] + 2.*x[0] - 2.
	gradient[1] = 200.*(x[1] - x[0]**2)
	return 0

def constraint_jacobian(number_variables, number_jacobian_nonzeros, x, jacobian_values, user_data):
	jacobian_values[0] = x[1]
	jacobian_values[1] = 1.
	jacobian_values[2] = x[0]
	jacobian_values[3] = 2.*x[1]
	return 0

def lagrangian_hessian(number_variables, number_constraints, number_hessian_nonzeros, x, objective_multiplier,
					   multipliers, hessian_values, user_data):
	hessian_values[0] = objective_multiplier*(1200*x[0]**2 - 400.*x[1] + 2.)
	hessian_values[1] = -400.*objective_multiplier*x[0] - multipliers[0]
	hessian_values[2] = 200.*objective_multiplier - 2.*multipliers[1]
	return 0

def lagrangian_hessian_operator(number_variables, number_constraints, x, evaluate_at_x, objective_multiplier, multipliers,
								vector, result, user_data):
	hessian00 = objective_multiplier*(1200*x[0]**2. - 400.*x[1] + 2.)
	hessian10 = -400.*objective_multiplier*x[0] - multipliers[0]
	hessian11 = 200.*objective_multiplier - 2.*multipliers[1]
	result[0] = hessian00*vector[0] + hessian10*vector[1]
	result[1] = hessian10*vector[0] + hessian11*vector[1]
	return 0

if __name__ == '__main__':
	# model creation
	problem_type = unopy.PROBLEM_NONLINEAR
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
	# user data
	user_data = {"offset": -306.5}

	model = unopy.Model(problem_type, number_variables, variables_lower_bounds, variables_upper_bounds, base_indexing)
	model.set_objective(optimization_sense, objective, objective_gradient)
	model.set_constraints(number_constraints, constraints, constraints_lower_bounds, constraints_upper_bounds,
	  number_jacobian_nonzeros, jacobian_row_indices, jacobian_column_indices, constraint_jacobian)
	model.set_lagrangian_hessian(number_hessian_nonzeros, hessian_triangular_part, hessian_row_indices,
		hessian_column_indices, lagrangian_hessian, lagrangian_sign_convention)
	#model.set_lagrangian_hessian_operator(number_hessian_nonzeros, lagrangian_hessian_operator, lagrangian_sign_convention)
	model.set_initial_primal_iterate(x0)
	model.set_user_data(user_data)
	
	uno_solver = unopy.UnoSolver()
	uno_solver.set_preset("filtersqp")
	#uno_solver.set_option("print_subproblem", "yes")
	#uno_solver.set_option("print_solution", "yes")
	#uno_solver.set_option("logger", "DEBUG3")
	
	result = uno_solver.optimize(model)

	# optimization summary
	print("\nReading optimization summary from Python:")
	print("Number of iterations:", result.optimization_status)
	print("Number of iterations:", result.solution_status)
	print("Objective at solution:", result.solution_objective)
	print("Primal feasibility at solution:", result.solution_primal_feasibility)
	print("Dual feasibility at solution:", result.solution_dual_feasibility)
	print("Complementarity at solution:", result.solution_complementarity)
	print("Primal solution:", result.primal_solution)
	print("Constraint dual solution:", result.constraint_dual_solution)
	print("Lower bound dual solution:", result.lower_bound_dual_solution)
	print("Upper bound dual solution:", result.upper_bound_dual_solution)
	print("Number of iterations:", result.number_iterations)
	print("CPU time:", result.cpu_time)
	print("Number of objective evaluations:", result.number_objective_evaluations)
	print("Number of constraint evaluations:", result.number_constraint_evaluations)
	print("Number of objective gradient evaluations:", result.number_objective_gradient_evaluations)
	print("Number of Jacobian evaluations:", result.number_jacobian_evaluations)
	print("Number of Hessian evaluations:", result.number_hessian_evaluations)
	print("Number of subproblems solved:", result.number_subproblems_solved)