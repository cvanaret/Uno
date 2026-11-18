// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Tests of the functions of the C API that get the result of an optimization and of the model functions, compiled as
// C: a C caller cannot catch C++ exceptions, so any exception that crossed the C interface would abort this program

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "uno/Uno_C_API.h"

static int number_failures = 0;

#define CHECK(condition) \
   do { \
      if (!(condition)) { \
         fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
         ++number_failures; \
      } \
   } while (0)

// min (x0 - 1)^2 + (x1 - 2)^2
static uno_int objective(uno_int number_variables, const double* x, double* objective_value, void* user_data) {
   (void)number_variables;
   (void)user_data;
   *objective_value = pow(x[0] - 1., 2.) + pow(x[1] - 2., 2.);
   return 0;
}

static uno_int objective_gradient(uno_int number_variables, const double* x, double* gradient, void* user_data) {
   (void)number_variables;
   (void)user_data;
   gradient[0] = 2. * (x[0] - 1.);
   gradient[1] = 2. * (x[1] - 2.);
   return 0;
}

static uno_int lagrangian_hessian(uno_int number_variables, uno_int number_constraints, uno_int number_hessian_nonzeros,
      const double* x, double objective_multiplier, const double* multipliers, double* hessian_values, void* user_data) {
   (void)number_variables;
   (void)number_constraints;
   (void)number_hessian_nonzeros;
   (void)x;
   (void)multipliers;
   (void)user_data;
   hessian_values[0] = 2. * objective_multiplier;
   hessian_values[1] = 2. * objective_multiplier;
   return 0;
}

static void check_not_available(void* solver) {
   uno_int integer_value = -1;
   double double_value = -1.;
   double vector[2] = {-1., -1.};
   CHECK(uno_get_optimization_status(solver, &integer_value) == UNO_RESULT_NOT_AVAILABLE);
   CHECK(uno_get_solution_objective(solver, &double_value) == UNO_RESULT_NOT_AVAILABLE);
   CHECK(uno_get_primal_solution_component(solver, 0, &double_value) == UNO_RESULT_NOT_AVAILABLE);
   CHECK(uno_get_primal_solution(solver, vector) == UNO_RESULT_NOT_AVAILABLE);
   CHECK(uno_get_number_iterations(solver, &integer_value) == UNO_RESULT_NOT_AVAILABLE);
   CHECK(uno_get_cpu_time(solver, &double_value) == UNO_RESULT_NOT_AVAILABLE);
   const char* description = NULL;
   CHECK(uno_get_method_description(solver, &description) == UNO_RESULT_NOT_AVAILABLE && description == NULL);
   CHECK(integer_value == -1 && double_value == -1. && vector[0] == -1. && vector[1] == -1.);
}

int main(void) {
   // invalid models
   CHECK(uno_create_unconstrained_model(NULL, 2, UNO_ZERO_BASED_INDEXING) == NULL);
   CHECK(uno_create_model(UNO_PROBLEM_NONLINEAR, 0, NULL, NULL, UNO_ZERO_BASED_INDEXING) == NULL);

   void* model = uno_create_unconstrained_model(UNO_PROBLEM_NONLINEAR, 2, UNO_ZERO_BASED_INDEXING);
   void* solver = uno_create_solver();
   CHECK(model != NULL && solver != NULL);
   const uno_int hessian_row_indices[] = {0, 1};
   const uno_int hessian_column_indices[] = {0, 1};
   CHECK(uno_set_objective(model, UNO_MINIMIZE, objective, objective_gradient));
   CHECK(!uno_set_lagrangian_hessian(model, 2, UNO_LOWER_TRIANGLE, NULL, NULL, lagrangian_hessian));
   CHECK(uno_set_lagrangian_hessian(model, 2, UNO_LOWER_TRIANGLE, hessian_row_indices, hessian_column_indices,
      lagrangian_hessian));
   CHECK(uno_set_solver_string_option(solver, "logger", "SILENT") == UNO_OPTIONS_SUCCESS);
   CHECK(uno_set_solver_preset(solver, "ipopt") == UNO_OPTIONS_SUCCESS);

   // no result before the solve
   check_not_available(solver);

   // null arguments
   double value = -1.;
   CHECK(uno_get_solution_objective(NULL, &value) == UNO_RESULT_INVALID_ARGUMENT);
   CHECK(uno_get_solution_objective(solver, NULL) == UNO_RESULT_INVALID_ARGUMENT);
   const char* description = NULL;
   CHECK(uno_get_method_description(NULL, &description) == UNO_RESULT_INVALID_ARGUMENT);

   // successful solve
   CHECK(uno_optimize(solver, model));
   uno_int optimization_status = -1;
   CHECK(uno_get_optimization_status(solver, &optimization_status) == UNO_RESULT_SUCCESS);
   CHECK(optimization_status == UNO_SUCCESS);
   double primal_solution[2];
   CHECK(uno_get_primal_solution(solver, primal_solution) == UNO_RESULT_SUCCESS);
   CHECK(fabs(primal_solution[0] - 1.) < 1e-6 && fabs(primal_solution[1] - 2.) < 1e-6);
   CHECK(uno_get_primal_solution_component(solver, 1, &value) == UNO_RESULT_SUCCESS);
   CHECK(value == primal_solution[1]);
   CHECK(uno_get_primal_solution_component(solver, 2, &value) == UNO_RESULT_INVALID_INDEX);
   CHECK(uno_get_primal_solution_component(solver, -1, &value) == UNO_RESULT_INVALID_INDEX);
   // no constraints: an empty vector can be retrieved with a null pointer
   CHECK(uno_get_constraint_dual_solution(solver, NULL) == UNO_RESULT_SUCCESS);
   CHECK(uno_get_constraint_dual_solution_component(solver, 0, &value) == UNO_RESULT_INVALID_INDEX);
   CHECK(uno_get_method_description(solver, &description) == UNO_RESULT_SUCCESS && description != NULL);

   // a failed optimization discards the result
   CHECK(!uno_optimize(solver, NULL));
   check_not_available(solver);

   uno_destroy_solver(solver);
   uno_destroy_model(model);
   if (number_failures > 0) {
      fprintf(stderr, "%d check(s) failed\n", number_failures);
      return EXIT_FAILURE;
   }
   printf("All result tests passed.\n");
   return EXIT_SUCCESS;
}
