#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "Uno_C_API.h"

int32_t objective_function(int32_t /*number_variables*/, const double* x, double* objective_value, void* /*user_data*/) {
	*objective_value = 100.*pow(x[1] - pow(x[0], 2.), 2.) + pow(1. - x[0], 2.);
	return 0;
}

int32_t constraint_functions(int32_t /*number_variables*/, int32_t /*number_constraints*/, const double* x,
		double* constraint_values, void* /*user_data*/) {
	constraint_values[0] = x[0]*x[1];
	constraint_values[1] = x[0] + pow(x[1], 2.);
	return 0;
}

int32_t objective_gradient(int32_t /*number_variables*/, const double* x, double* gradient, void* /*user_data*/) {
	gradient[0] = 400.*pow(x[0], 3.) - 400.*x[0]*x[1] + 2.*x[0] - 2.;
	gradient[1] = 200.*(x[1] - pow(x[0], 2.));
	return 0;
}

int32_t constraint_jacobian(int32_t /*number_variables*/, int32_t /*number_jacobian_nonzeros*/, const double* x,
		double* jacobian, void* /*user_data*/) {
	jacobian[0] = x[1];
	jacobian[1] = 1.;
	jacobian[2] = x[0];
	jacobian[3] = 2.*x[1];
	return 0.;
}

int32_t lagrangian_hessian(int32_t /*number_variables*/, int32_t /*number_constraints*/, int32_t /*number_hessian_nonzeros*/,
      const double* x, double objective_multiplier, const double* multipliers, double* hessian, void* /*user_data*/) {
	hessian[0] = objective_multiplier*(1200*pow(x[0], 2.) - 400.*x[1] + 2.);
	hessian[1] = -400.*objective_multiplier*x[0] - multipliers[0];
	hessian[2] = 200.*objective_multiplier - 2.*multipliers[1];
	return 0;
}

int32_t lagrangian_hessian_operator(int32_t number_variables, int32_t number_constraints, const double* x,
      bool evaluate_at_x, double objective_multiplier, const double* multipliers, const double* vector,
      double* result, void* user_data) {
	const double hessian00 = objective_multiplier*(1200*pow(x[0], 2.) - 400.*x[1] + 2.);
	const double hessian10 = -400.*objective_multiplier*x[0] - multipliers[0];
	const double hessian11 = 200.*objective_multiplier - 2.*multipliers[1];
	result[0] = hessian00*vector[0] + hessian10*vector[1];
	result[1] = hessian10*vector[0] + hessian11*vector[1];
	return 0;
}

void print_vector(const double* vector, int32_t size) {
	for (size_t index = 0; index < size; ++index) {
		printf("%g ", vector[index]);
	}
	printf("\n");
}

int main() {
	int32_t uno_major, uno_minor, uno_patch;
	uno_get_version(&uno_major, &uno_minor, &uno_patch);
	printf("Uno v%d.%d.%d\n", uno_major, uno_minor, uno_patch);

	// model creation
	const int32_t base_indexing = UNO_ZERO_BASED_INDEXING;
	// variables
	const int32_t number_variables = 2;
	double variables_lower_bounds[] = {-INFINITY, -INFINITY};
	double variables_upper_bounds[] = {0.5, INFINITY};
	// objective
	const int32_t optimization_sense = UNO_MINIMIZE;
	// constraints
	const int32_t number_constraints = 2;
	double constraints_lower_bounds[] = {1., 0.};
	double constraints_upper_bounds[] = {INFINITY, INFINITY};
	const int32_t number_jacobian_nonzeros = 4;
	int32_t jacobian_row_indices[] = {0, 1, 0, 1};
	int32_t jacobian_column_indices[] = {0, 0, 1, 1};
	// Hessian
	const int32_t number_hessian_nonzeros = 3;
	const char hessian_triangular_part = UNO_LOWER_TRIANGLE;
	int32_t hessian_row_indices[] = {0, 1, 1};
	int32_t hessian_column_indices[] = {0, 0, 1};
	const double lagrangian_sign_convention = UNO_MULTIPLIER_NEGATIVE;
	// initial point
	double x0[] = {-2., 1.};

	void* model = uno_create_model(UNO_PROBLEM_NONLINEAR, number_variables, variables_lower_bounds,
      variables_upper_bounds, base_indexing);
	assert(uno_set_objective(model, optimization_sense, objective_function, objective_gradient));
	assert(uno_set_constraints(model, number_constraints, constraint_functions,
		constraints_lower_bounds, constraints_upper_bounds, number_jacobian_nonzeros,
		jacobian_row_indices, jacobian_column_indices, constraint_jacobian));
	assert(uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part, hessian_row_indices,
		hessian_column_indices, lagrangian_hessian, lagrangian_sign_convention));
/*
	assert(uno_set_lagrangian_hessian_operator(model, number_hessian_nonzeros, lagrangian_hessian_operator,
		lagrangian_sign_convention));
*/
	assert(uno_set_initial_primal_iterate(model, x0));

	// solver creation
	void* solver = uno_create_solver();
	uno_set_solver_preset(solver, "filtersqp");
	uno_set_solver_option(solver, "print_solution", "yes");

	// solve
	uno_optimize(solver, model);

	// get the solution
	const int32_t optimization_status = uno_get_optimization_status(solver);
	assert(optimization_status == UNO_SUCCESS);
	const int32_t iterate_status = uno_get_solution_status(solver);
	assert(iterate_status == UNO_FEASIBLE_KKT_POINT);
	const double solution_objective = uno_get_solution_objective(solver);
	printf("\nReading the result from the C file:\n");
	printf("Solution objective = %g\n", solution_objective);

	double primal_solution[number_variables];
	uno_get_primal_solution(solver, primal_solution);
	printf("Primal solution: "); print_vector(primal_solution, number_variables);
	double constraint_dual_solution[number_constraints];
	uno_get_constraint_dual_solution(solver, constraint_dual_solution);
	printf("Constraint dual solution: "); print_vector(constraint_dual_solution, number_constraints);
	double lower_bound_dual_solution[number_variables];
	uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution);
	printf("Lower bound dual solution: "); print_vector(lower_bound_dual_solution, number_variables);
	double upper_bound_dual_solution[number_variables];
	uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution);
	printf("Upper bound dual solution: "); print_vector(upper_bound_dual_solution, number_variables);
	const double solution_primal_feasibility = uno_get_solution_primal_feasibility(solver);
	printf("Primal feasibility solution at solution = %e\n", solution_primal_feasibility);
	const double solution_dual_feasibility = uno_get_solution_dual_feasibility(solver);
	printf("Dual feasibility solution at solution = %e\n", solution_dual_feasibility);
	const double solution_complementarity = uno_get_solution_complementarity(solver);
	printf("Complementarity solution at solution = %e\n", solution_complementarity);

	// cleanup
	uno_destroy_solver(solver);
	uno_destroy_model(model);
	return 0;
}