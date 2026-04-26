# Copyright (c) 2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

import unopy
Inf = float("inf")

# hs015.mod
def objective(x):
	# the modeler should throw an exception if the function cannot be evaluated
	return 100.*(x[1] - x[0]**2)**2 + (1. - x[0])**2
	
def constraints(x, constraint_values):
	# the modeler should throw an exception if the function cannot be evaluated
	constraint_values[:] = [x[0]*x[1], x[0] + x[1]**2]
	
def objective_gradient(x, gradient):
	# the modeler should throw an exception if the function cannot be evaluated
	gradient[:] = [400.*x[0]**3 - 400.*x[0]*x[1] + 2.*x[0] - 2.,
				   200.*(x[1] - x[0]**2)]

def jacobian(x, jacobian_values):
	# the modeler should throw an exception if the function cannot be evaluated
	jacobian_values[:] = [x[1], 1., x[0], 2.*x[1]]

def lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values):
	# the modeler should throw an exception if the function cannot be evaluated
	hessian_values[:] = [objective_multiplier*(1200*x[0]**2 - 400.*x[1] + 2.),
						-400.*objective_multiplier*x[0] - multipliers[0],
						200.*objective_multiplier - 2.*multipliers[1]]

def lagrangian_hessian_operator(x, evaluate_at_x, objective_multiplier, multipliers, vector, result):
	# the modeler should throw an exception if the function cannot be evaluated
	hessian00 = objective_multiplier*(1200*x[0]**2. - 400.*x[1] + 2.)
	hessian10 = -400.*objective_multiplier*x[0] - multipliers[0]
	hessian11 = 200.*objective_multiplier - 2.*multipliers[1]
	result[:] = [hessian00*vector[0] + hessian10*vector[1],
				hessian10*vector[0] + hessian11*vector[1]]

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

	model = unopy.Model(problem_type, number_variables, base_indexing)
	model.set_variables_lower_bounds(variables_lower_bounds)
	model.set_variables_upper_bounds(variables_upper_bounds)
	model.set_objective(optimization_sense, objective, objective_gradient)
	model.set_constraints(number_constraints, constraints, constraints_lower_bounds, constraints_upper_bounds,
	  number_jacobian_nonzeros, jacobian_row_indices, jacobian_column_indices, jacobian)
	model.set_initial_primal_iterate(x0)

	# solver creation
	uno_solver = unopy.UnoSolver()
	uno_solver.set_preset("filtersqp")
	uno_solver.set_option("QP_solver", "BQPD")
	print("Solving with Uno", unopy.current_uno_version())

	# run 1: solve with the filtersqp preset with no exact Hessian. Uno defaults to L-BFGS Hessian for NLPs
	result = uno_solver.optimize(model)
	print("Objective at solution:", result.solution_objective)

	# run 2: solve with the filtersqp preset with exact Hessian
	model.set_lagrangian_hessian(number_hessian_nonzeros, hessian_triangular_part, hessian_row_indices,
	  hessian_column_indices, lagrangian_hessian)
	model.set_lagrangian_sign_convention(lagrangian_sign_convention)
	# the Hessian model was overwritten. Set it again
	uno_solver.set_option("hessian_model", "exact")
	result = uno_solver.optimize(model)
	print("Objective at solution:", result.solution_objective)
	# optimization summary
	print("Optimization status:", result.optimization_status)
	print("Solution status:", result.solution_status)
	print("Primal feasibility at solution:", result.solution_primal_feasibility)
	print("Stationarity at solution:", result.solution_stationarity)
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
	assert abs(result.solution_objective - 306.5) <= 1e-4

	# run 3: solve with the ipopt preset
	uno_solver.set_preset("ipopt")
	uno_solver.set_option("linear_solver", "MUMPS")
	result = uno_solver.optimize(model)
	assert abs(result.solution_objective - 306.5) <= 1e-4
	
	# run 4: solve with the filterslp preset
	uno_solver.set_preset("filterslp")
	uno_solver.set_option("LP_solver", "HiGHS")
	result = uno_solver.optimize(model)
	assert abs(result.solution_objective - 306.5) <= 1e-4