#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "uno/Uno_C_API.h"

uno_int objective_function(uno_int /*number_variables*/, const double* x, double* objective_value, void* /*user_data*/) {
   *objective_value = 100.*pow(x[1] - pow(x[0], 2.), 2.) + pow(1. - x[0], 2.);
   return 0;
}

uno_int constraint_functions(uno_int /*number_variables*/, uno_int /*number_constraints*/, const double* x,
      double* constraint_values, void* /*user_data*/) {
   constraint_values[0] = x[0]*x[1];
   constraint_values[1] = x[0] + pow(x[1], 2.);
   return 0;
}

uno_int objective_gradient(uno_int /*number_variables*/, const double* x, double* gradient, void* /*user_data*/) {
   gradient[0] = 400.*pow(x[0], 3.) - 400.*x[0]*x[1] + 2.*x[0] - 2.;
   gradient[1] = 200.*(x[1] - pow(x[0], 2.));
   return 0;
}

uno_int jacobian(uno_int /*number_variables*/, uno_int /*number_jacobian_nonzeros*/, const double* x,
      double* jacobian_values, void* /*user_data*/) {
   jacobian_values[0] = x[1];
   jacobian_values[1] = 1.;
   jacobian_values[2] = x[0];
   jacobian_values[3] = 2.*x[1];
   return 0;
}

uno_int jacobian_operator(uno_int /*number_variables*/, uno_int /*number_constraints*/, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* /*user_data*/) {
   result[0] = x[1]*vector[0] + 1.*vector[1];
   result[1] = x[0]*vector[0] + 2.*x[1]*vector[1];
   return 0;
}

uno_int jacobian_transposed_operator(uno_int /*number_variables*/, uno_int /*number_constraints*/, const double* x,
      bool evaluate_at_x, const double* vector, double* result, void* /*user_data*/) {
   result[0] = x[1]*vector[0] + x[0]*vector[1];
   result[1] = 1.*vector[0] + 2.*x[1]*vector[1];
   return 0;
}

uno_int lagrangian_hessian(uno_int /*number_variables*/, uno_int /*number_constraints*/, uno_int /*number_hessian_nonzeros*/,
      const double* x, double objective_multiplier, const double* multipliers, double* hessian_values, void* /*user_data*/) {
   hessian_values[0] = objective_multiplier*(1200*pow(x[0], 2.) - 400.*x[1] + 2.);
   hessian_values[1] = -400.*objective_multiplier*x[0] - multipliers[0];
   hessian_values[2] = 200.*objective_multiplier - 2.*multipliers[1];
   return 0;
}

uno_int lagrangian_hessian_operator(uno_int number_variables, uno_int number_constraints, const double* x,
      bool evaluate_at_x, double objective_multiplier, const double* multipliers, const double* vector,
      double* result, void* user_data) {
   const double hessian00 = objective_multiplier*(1200*pow(x[0], 2.) - 400.*x[1] + 2.);
   const double hessian10 = -400.*objective_multiplier*x[0] - multipliers[0];
   const double hessian11 = 200.*objective_multiplier - 2.*multipliers[1];
   result[0] = hessian00*vector[0] + hessian10*vector[1];
   result[1] = hessian10*vector[0] + hessian11*vector[1];
   return 0;
}

// stops if the result of the optimization could not be retrieved
void check_result(uno_int status, const char* function_name) {
   if (status != UNO_RESULT_SUCCESS) {
      fprintf(stderr, "%s failed with status %d\n", function_name, (int)status);
      exit(EXIT_FAILURE);
   }
}

const char* method_description(void* solver) {
   const char* description = NULL;
   check_result(uno_get_method_description(solver, &description), "uno_get_method_description");
   return description;
}

void print_vector(const double* vector, uno_int size) {
   for (size_t index = 0; index < size; ++index) {
      printf("%g ", vector[index]);
   }
   printf("\n");
}

void print_uno_options() {
   // temporary solver for default options
   void* solver = uno_create_solver();
   printf("Uno options:\n");
   for (uno_options_iterator it = uno_options_begin_iterator(); it != uno_options_end_iterator(); uno_options_next_iterator(&it)) {
      const char* name = uno_options_iterator_name(it);
      uno_int type = uno_options_iterator_type(it);
      uno_int status = UNO_OPTIONS_ERROR;
      if (type == UNO_OPTION_TYPE_DOUBLE) {
         double value;
         status = uno_get_solver_double_option(solver, name, &value);
         if (status == UNO_OPTIONS_SUCCESS) {
            printf("\t%s: %g\n", name, value);
         }
      }
      else if (type == UNO_OPTION_TYPE_INTEGER) {
         uno_int value;
         status = uno_get_solver_integer_option(solver, name, &value);
         if (status == UNO_OPTIONS_SUCCESS) {
            printf("\t%s: %d\n", name, value);
         }
      }
      else if (type == UNO_OPTION_TYPE_STRING) {
         const char* value;
         status = uno_get_solver_string_option(solver, name, &value);
         if (status == UNO_OPTIONS_SUCCESS) {
            printf("\t%s: %s\n", name, value);
         }
      }
      else if (type == UNO_OPTION_TYPE_BOOL) {
         bool value;
         status = uno_get_solver_bool_option(solver, name, &value);
         if (status == UNO_OPTIONS_SUCCESS) {
            printf("\t%s: %s\n", name, value ? "true" : "false");
         }
      }
      // some options have no value before the solve (e.g. the options set by the "auto" preset)
      if (status == UNO_OPTIONS_NOT_SET) {
         printf("\t%s: (not set)\n", name);
      }
      assert(status == UNO_OPTIONS_SUCCESS || status == UNO_OPTIONS_NOT_SET);
   }
   uno_destroy_solver(solver);
}

int main() {
   uno_int uno_major, uno_minor, uno_patch;
   uno_get_version(&uno_major, &uno_minor, &uno_patch);
   printf("Uno v%d.%d.%d\n", uno_major, uno_minor, uno_patch);

   // print all options
   print_uno_options();

   // model creation
   const uno_int base_indexing = UNO_ZERO_BASED_INDEXING;
   // variables
   const uno_int number_variables = 2;
   double variables_lower_bounds[] = {-INFINITY, -INFINITY};
   double variables_upper_bounds[] = {0.5, INFINITY};
   // objective
   const uno_int optimization_sense = UNO_MINIMIZE;
   // constraints
   const uno_int number_constraints = 2;
   double constraints_lower_bounds[] = {1., 0.};
   double constraints_upper_bounds[] = {INFINITY, INFINITY};
   const uno_int number_jacobian_nonzeros = 4;
   uno_int jacobian_row_indices[] = {0, 1, 0, 1};
   uno_int jacobian_column_indices[] = {0, 0, 1, 1};
   // Hessian
   const uno_int number_hessian_nonzeros = 3;
   const char hessian_triangular_part = UNO_LOWER_TRIANGLE;
   uno_int hessian_row_indices[] = {0, 1, 1};
   uno_int hessian_column_indices[] = {0, 0, 1};
   const uno_int lagrangian_sign_convention = UNO_MULTIPLIER_NEGATIVE;
   // initial point
   double x0[] = {-2., 1.};

   void* model = uno_create_model(UNO_PROBLEM_NONLINEAR, number_variables, variables_lower_bounds,
      variables_upper_bounds, base_indexing);
   uno_set_model_name(model, "hs015");
   uno_set_objective(model, optimization_sense, objective_function, objective_gradient);
   uno_set_constraints(model, number_constraints, constraint_functions,
      constraints_lower_bounds, constraints_upper_bounds, number_jacobian_nonzeros,
      jacobian_row_indices, jacobian_column_indices, jacobian);
/*
   uno_set_jacobian_operator(model, jacobian_operator));
   uno_set_jacobian_transposed_operator(model, jacobian_transposed_operator));
   uno_set_lagrangian_hessian_operator(model, lagrangian_hessian_operator));
*/
   uno_set_initial_primal_iterate(model, x0);

   // solver creation
   void* solver = uno_create_solver();
   uno_int option_status = uno_set_solver_preset(solver, "filtersqp");
   assert(option_status == UNO_OPTIONS_SUCCESS);
   option_status = uno_set_solver_bool_option(solver, "print_solution", true);
   assert(option_status == UNO_OPTIONS_SUCCESS);

   // run 1: solve with no Hessian. Uno defaults to L-BFGS Hessian for NLPs
   uno_int optimize_status = uno_optimize(solver, model);
   assert(optimize_status == UNO_OPTIMIZE_SUCCESS);
   printf("Method used: %s\n", method_description(solver));
   // get the solution
   uno_int optimization_status;
   check_result(uno_get_optimization_status(solver, &optimization_status), "uno_get_optimization_status");
   assert(optimization_status == UNO_SUCCESS);
   uno_int iterate_status;
   check_result(uno_get_solution_status(solver, &iterate_status), "uno_get_solution_status");
   assert(iterate_status == UNO_FEASIBLE_KKT_POINT);
   double solution_objective;
   check_result(uno_get_solution_objective(solver, &solution_objective), "uno_get_solution_objective");
   printf("Solution objective = %g\n", solution_objective);
   const char* hessian_model = NULL;
   option_status = uno_get_solver_string_option(solver, "hessian_model", &hessian_model);
   assert(option_status == UNO_OPTIONS_SUCCESS);
   if (hessian_model != NULL) {
      printf("Hessian model picked by Uno = %s\n", hessian_model);
      assert(strcmp(hessian_model, "LBFGS") == 0);
   }
   printf("\n");

   // run 2: solve with exact Hessian
   uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part, hessian_row_indices,
      hessian_column_indices, lagrangian_hessian);
   uno_set_lagrangian_sign_convention(model, lagrangian_sign_convention);
   optimize_status = uno_optimize(solver, model);
   assert(optimize_status == UNO_OPTIMIZE_SUCCESS);
   
   printf("Method used: %s\n", method_description(solver));
   // get the solution
   check_result(uno_get_optimization_status(solver, &optimization_status), "uno_get_optimization_status");
   assert(optimization_status == UNO_SUCCESS);
   check_result(uno_get_solution_status(solver, &iterate_status), "uno_get_solution_status");
   assert(iterate_status == UNO_FEASIBLE_KKT_POINT);
   check_result(uno_get_solution_objective(solver, &solution_objective), "uno_get_solution_objective");
   printf("Solution objective = %g\n", solution_objective);
   double primal_solution[number_variables];
   check_result(uno_get_primal_solution(solver, primal_solution), "uno_get_primal_solution");
   printf("Primal solution: "); print_vector(primal_solution, number_variables);
   double constraint_dual_solution[number_constraints];
   check_result(uno_get_constraint_dual_solution(solver, constraint_dual_solution), "uno_get_constraint_dual_solution");
   printf("Constraint dual solution: "); print_vector(constraint_dual_solution, number_constraints);
   double lower_bound_dual_solution[number_variables];
   check_result(uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution), "uno_get_lower_bound_dual_solution");
   printf("Lower bound dual solution: "); print_vector(lower_bound_dual_solution, number_variables);
   double upper_bound_dual_solution[number_variables];
   check_result(uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution), "uno_get_upper_bound_dual_solution");
   printf("Upper bound dual solution: "); print_vector(upper_bound_dual_solution, number_variables);
   double solution_constraints[number_constraints];
   check_result(uno_get_solution_constraints(solver, solution_constraints), "uno_get_solution_constraints");
   printf("Solution constraints: "); print_vector(solution_constraints, number_constraints);
   
   double solution_primal_feasibility;
   check_result(uno_get_solution_primal_feasibility(solver, &solution_primal_feasibility), "uno_get_solution_primal_feasibility");
   printf("Primal feasibility at solution = %e\n", solution_primal_feasibility);
   double solution_stationarity;
   check_result(uno_get_solution_stationarity(solver, &solution_stationarity), "uno_get_solution_stationarity");
   printf("Stationarity at solution = %e\n", solution_stationarity);
   double solution_complementarity;
   check_result(uno_get_solution_complementarity(solver, &solution_complementarity), "uno_get_solution_complementarity");
   printf("Complementarity at solution = %e\n", solution_complementarity);

#ifdef UNO_EXAMPLE_HIGHS
   // run 3: solve with the filterslp preset and HiGHS as LP solver
   printf("\n");
   option_status = uno_set_solver_preset(solver, "filterslp");
   assert(option_status == UNO_OPTIONS_SUCCESS);
   option_status = uno_set_solver_string_option(solver, "LP_solver", "HIGHS");
   assert(option_status == UNO_OPTIONS_SUCCESS);
   optimize_status = uno_optimize(solver, model);
   assert(optimize_status == UNO_OPTIMIZE_SUCCESS);
   printf("Method used: %s\n", method_description(solver));
   check_result(uno_get_optimization_status(solver, &optimization_status), "uno_get_optimization_status");
   assert(optimization_status == UNO_SUCCESS);
   check_result(uno_get_solution_status(solver, &iterate_status), "uno_get_solution_status");
   assert(iterate_status == UNO_FEASIBLE_KKT_POINT);
   check_result(uno_get_solution_objective(solver, &solution_objective), "uno_get_solution_objective");
   printf("Solution objective = %g\n", solution_objective);
#endif

   // cleanup
   uno_destroy_solver(solver);
   uno_destroy_model(model);
   return 0;
}