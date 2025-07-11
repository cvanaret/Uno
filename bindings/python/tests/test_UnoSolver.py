# Copyright (c) 2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

import unopy

def objective(x):
	return x[0] + x[1]

if __name__ == '__main__':
	options = unopy.Options.get_default()
	unopy.Options.set_preset(options, "filtersqp")
	solver = unopy.UnoSolver(True, options)
	number_variables = 2
	number_constraints = 0
	solver.solve(number_variables, number_constraints, objective)
