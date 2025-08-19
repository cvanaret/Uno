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
	uno_set_objective(model, optimization_sense, objective_function, objective_gradient);
	uno_set_constraints(model, number_constraints, constraint_functions,
		constraints_lower_bounds, constraints_upper_bounds, number_jacobian_nonzeros,
		jacobian_row_indices, jacobian_column_indices, constraint_jacobian);
	uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part, hessian_row_indices,
		hessian_column_indices, lagrangian_hessian, lagrangian_sign_convention);
	uno_set_initial_primal_iterate(model, x0);

	// solver creation
	void* solver = uno_create_solver();
	// uno_set_solver_preset(solver, "ipopt");

	// solve
	uno_optimize(solver, model);

	// cleanup
	uno_destroy_solver(solver);
	uno_destroy_model(model);
	return 0;
}