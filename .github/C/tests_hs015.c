#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "Uno_C_API.h"

int32_t objective_function(int32_t /*number_variables*/, const double* x, double* objective_value, void* /*user_data*/) {
	*objective_value = 100.*pow(x[1] - pow(x[0], 2.), 2.) + pow(1. - x[0], 2.);
	return 0;
}

int32_t objective_function_max(int32_t /*number_variables*/, const double* x, double* objective_value, void* /*user_data*/) {
	*objective_value = -(100.*pow(x[1] - pow(x[0], 2.), 2.) + pow(1. - x[0], 2.));
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

int32_t objective_gradient_max(int32_t /*number_variables*/, const double* x, double* gradient, void* /*user_data*/) {
	gradient[0] = -(400.*pow(x[0], 3.) - 400.*x[0]*x[1] + 2.*x[0] - 2.);
	gradient[1] = -(200.*(x[1] - pow(x[0], 2.)));
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

int32_t lagrangian_hessian_negative_sign(int32_t /*number_variables*/, int32_t /*number_constraints*/, int32_t /*number_hessian_nonzeros*/,
      const double* x, double objective_multiplier, const double* multipliers, double* hessian, void* /*user_data*/) {
	hessian[0] = objective_multiplier*(1200*pow(x[0], 2.) - 400.*x[1] + 2.);
	hessian[1] = -400.*objective_multiplier*x[0] - multipliers[0];
	hessian[2] = 200.*objective_multiplier - 2.*multipliers[1];
	return 0;
}

int32_t lagrangian_hessian_positive_sign(int32_t /*number_variables*/, int32_t /*number_constraints*/, int32_t /*number_hessian_nonzeros*/,
      const double* x, double objective_multiplier, const double* multipliers, double* hessian, void* /*user_data*/) {
	hessian[0] = objective_multiplier*(1200*pow(x[0], 2.) - 400.*x[1] + 2.);
	hessian[1] = -400.*objective_multiplier*x[0] + multipliers[0];
	hessian[2] = 200.*objective_multiplier + 2.*multipliers[1];
	return 0;
}

int32_t lagrangian_hessian_max_negative_sign(int32_t /*number_variables*/, int32_t /*number_constraints*/, int32_t /*number_hessian_nonzeros*/,
		const double* x, double objective_multiplier, const double* multipliers, double* hessian, void* /*user_data*/) {
	hessian[0] = -objective_multiplier*(1200*pow(x[0], 2.) - 400.*x[1] + 2.);
	hessian[1] = -400.*-objective_multiplier*x[0] - multipliers[0];
	hessian[2] = 200.*-objective_multiplier - 2.*multipliers[1];
	return 0;
}

int32_t lagrangian_hessian_max_positive_sign(int32_t /*number_variables*/, int32_t /*number_constraints*/, int32_t /*number_hessian_nonzeros*/,
		const double* x, double objective_multiplier, const double* multipliers, double* hessian, void* /*user_data*/) {
	hessian[0] = -objective_multiplier*(1200*pow(x[0], 2.) - 400.*x[1] + 2.);
	hessian[1] = -400.*-objective_multiplier*x[0] + multipliers[0];
	hessian[2] = 200.*-objective_multiplier + 2.*multipliers[1];
	return 0;
}

void print_vector(const double* vector, int32_t size) {
	for (size_t index = 0; index < size; ++index) {
		printf("%g ", vector[index]);
	}
	printf("\n");
}

void solve_instance(int32_t optimization_sense, double lagrangian_sign_convention, double reference_objective,
		const double reference_primal_solution[], const double reference_constraint_dual_solution[],
		const double reference_lower_bound_dual_solution[], const double reference_upper_bound_dual_solution[]) {
	// model creation
	const int32_t base_indexing = UNO_ZERO_BASED_INDEXING;
	// variables
	const int32_t number_variables = 2;
	double variables_lower_bounds[] = {-INFINITY, -INFINITY};
	double variables_upper_bounds[] = {0.5, INFINITY};
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
	// initial point
	double x0[] = {-2., 1.};
	
	void* model = uno_create_model(UNO_PROBLEM_NONLINEAR, number_variables, variables_lower_bounds,
      variables_upper_bounds, base_indexing);
	if (optimization_sense == UNO_MINIMIZE) {
		assert(uno_set_objective(model, optimization_sense, objective_function, objective_gradient));
	}
	else {
		assert(uno_set_objective(model, optimization_sense, objective_function_max, objective_gradient_max));
	}
	assert(uno_set_constraints(model, number_constraints, constraint_functions,
		constraints_lower_bounds, constraints_upper_bounds, number_jacobian_nonzeros,
		jacobian_row_indices, jacobian_column_indices, constraint_jacobian));
	if (lagrangian_sign_convention == UNO_MULTIPLIER_NEGATIVE) {
		if (optimization_sense == UNO_MINIMIZE) {
			assert(uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part,
				hessian_row_indices, hessian_column_indices, lagrangian_hessian_negative_sign,
				lagrangian_sign_convention));
		}
		else {
			assert(uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part,
				hessian_row_indices, hessian_column_indices, lagrangian_hessian_max_negative_sign,
				lagrangian_sign_convention));
		}
	}
	else {
		if (optimization_sense == UNO_MINIMIZE) {
			assert(uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part,
				hessian_row_indices, hessian_column_indices, lagrangian_hessian_positive_sign,
				lagrangian_sign_convention));
		}
		else {
			assert(uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part,
				hessian_row_indices, hessian_column_indices, lagrangian_hessian_max_positive_sign,
				lagrangian_sign_convention));
		}
	}
	assert(uno_set_initial_primal_iterate(model, x0));

	// solver creation
	void* solver = uno_create_solver();
	uno_set_solver_preset(solver, "filtersqp");
	uno_set_solver_option(solver, "logger", "SILENT");

	// solve
	uno_optimize(solver, model);
	
	// tests
	const int32_t optimization_status = uno_get_optimization_status(solver);
	assert(optimization_status == UNO_SUCCESS);
	const int32_t iterate_status = uno_get_solution_status(solver);
	assert(iterate_status == UNO_FEASIBLE_KKT_POINT);
	// objective
	const double solution_objective = uno_get_solution_objective(solver);
	assert(solution_objective == reference_objective);
	// primal solution
	double primal_solution[number_variables];
	uno_get_primal_solution(solver, primal_solution);
	for (size_t variable_index = 0; variable_index < number_variables; ++variable_index) {
		assert(primal_solution[variable_index] == reference_primal_solution[variable_index]);
	}
	// constraint dual solution
	double constraint_dual_solution[number_constraints];
	uno_get_constraint_dual_solution(solver, constraint_dual_solution);
	for (size_t constraint_index = 0; constraint_index < number_constraints; ++constraint_index) {
		assert(constraint_dual_solution[constraint_index] == reference_constraint_dual_solution[constraint_index]);
	}
	// lower bound dual solution
	double lower_bound_dual_solution[number_variables];
	uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution);
	for (size_t variable_index = 0; variable_index < number_variables; ++variable_index) {
		assert(lower_bound_dual_solution[variable_index] == reference_lower_bound_dual_solution[variable_index]);
	}
	// upper bound dual solution
	double upper_bound_dual_solution[number_variables];
	uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution);
	for (size_t variable_index = 0; variable_index < number_variables; ++variable_index) {
		assert(upper_bound_dual_solution[variable_index] == reference_upper_bound_dual_solution[variable_index]);
	}
	
	// cleanup
	uno_destroy_solver(solver);
	uno_destroy_model(model);
}

int main() {
	const double reference_objective = 306.5;
	const double reference_primal_solution[] = {0.5, 2};
	const double reference_constraint_dual_solution_negative[] = {700., 0.};
	const double reference_constraint_dual_solution_positive[] = {-700., 0.};
	const double reference_lower_bound_dual_solution[] = {0., 0.};
	const double reference_upper_bound_dual_solution_negative[] = {-1751., 0.};
	const double reference_upper_bound_dual_solution_positive[] = {1751., 0.};
	
	solve_instance(UNO_MINIMIZE, UNO_MULTIPLIER_NEGATIVE, reference_objective,
		reference_primal_solution, reference_constraint_dual_solution_negative,
		reference_lower_bound_dual_solution, reference_upper_bound_dual_solution_negative);
	printf("(UNO_MINIMIZE, UNO_MULTIPLIER_NEGATIVE) passed.\n");
	solve_instance(UNO_MINIMIZE, UNO_MULTIPLIER_POSITIVE, reference_objective,
		reference_primal_solution, reference_constraint_dual_solution_positive,
		reference_lower_bound_dual_solution, reference_upper_bound_dual_solution_positive);
	printf("(UNO_MINIMIZE, UNO_MULTIPLIER_POSITIVE) passed.\n");
	solve_instance(UNO_MAXIMIZE, UNO_MULTIPLIER_NEGATIVE, -reference_objective,
			reference_primal_solution, reference_constraint_dual_solution_positive,
			reference_lower_bound_dual_solution, reference_upper_bound_dual_solution_positive);
	printf("(UNO_MAXIMIZE, UNO_MULTIPLIER_NEGATIVE) passed.\n");
	solve_instance(UNO_MAXIMIZE, UNO_MULTIPLIER_POSITIVE, -reference_objective,
		reference_primal_solution, reference_constraint_dual_solution_negative,
		reference_lower_bound_dual_solution, reference_upper_bound_dual_solution_negative);
	printf("(UNO_MAXIMIZE, UNO_MULTIPLIER_POSITIVE) passed.\n");
	
	return 0;
}