// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Tests of the functions of the C API that deal with options, compiled as C: a C caller cannot catch C++ exceptions,
// so any exception that crossed the C interface would abort this program

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uno/Uno_C_API.h"

static int number_failures = 0;

#define CHECK(condition) \
   do { \
      if (!(condition)) { \
         fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
         ++number_failures; \
      } \
   } while (0)

static void test_setters(void) {
   void* solver = uno_create_solver();
   CHECK(uno_set_solver_integer_option(solver, "max_iterations", 50) == UNO_OPTIONS_SUCCESS);
   CHECK(uno_set_solver_double_option(solver, "primal_tolerance", 1e-6) == UNO_OPTIONS_SUCCESS);
   CHECK(uno_set_solver_bool_option(solver, "print_solution", true) == UNO_OPTIONS_SUCCESS);
   CHECK(uno_set_solver_string_option(solver, "hessian_model", "LBFGS") == UNO_OPTIONS_SUCCESS);
   CHECK(uno_set_solver_preset(solver, "filtersqp") == UNO_OPTIONS_SUCCESS);

   // unknown option
   CHECK(uno_set_solver_integer_option(solver, "max_iteration", 50) == UNO_OPTIONS_NOT_FOUND);
   CHECK(uno_set_solver_string_option(solver, "loger", "INFO") == UNO_OPTIONS_NOT_FOUND);
   // wrong type
   CHECK(uno_set_solver_double_option(solver, "max_iterations", 50.) == UNO_OPTIONS_WRONG_TYPE);
   CHECK(uno_set_solver_bool_option(solver, "logger", true) == UNO_OPTIONS_WRONG_TYPE);
   // invalid value (#1070)
   CHECK(uno_set_solver_string_option(solver, "logger", "SOLENT") == UNO_OPTIONS_INVALID_VALUE);
   CHECK(uno_set_solver_preset(solver, "snopt") == UNO_OPTIONS_INVALID_VALUE);
   // null arguments
   CHECK(uno_set_solver_integer_option(NULL, "max_iterations", 50) == UNO_OPTIONS_INVALID_ARGUMENT);
   CHECK(uno_set_solver_integer_option(solver, NULL, 50) == UNO_OPTIONS_INVALID_ARGUMENT);
   CHECK(uno_set_solver_string_option(solver, "logger", NULL) == UNO_OPTIONS_INVALID_ARGUMENT);
   CHECK(uno_set_solver_preset(solver, NULL) == UNO_OPTIONS_INVALID_ARGUMENT);
   // option file
   CHECK(uno_load_solver_option_file(solver, "this_option_file_does_not_exist.opt") == UNO_OPTIONS_INVALID_FILE);
   CHECK(uno_load_solver_option_file(solver, NULL) == UNO_OPTIONS_INVALID_ARGUMENT);
   // option type
   CHECK(uno_get_solver_option_type(solver, "logger") == UNO_OPTION_TYPE_STRING);
   CHECK(uno_get_solver_option_type(solver, "loger") == UNO_OPTION_TYPE_NOT_FOUND);
   CHECK(uno_get_solver_option_type(NULL, "logger") == UNO_OPTION_TYPE_NOT_FOUND);
   uno_destroy_solver(solver);
}

static void test_optimize_with_invalid_arguments(void) {
   void* solver = uno_create_solver();
   void* model = uno_create_unconstrained_model(UNO_PROBLEM_NONLINEAR, 2, UNO_ZERO_BASED_INDEXING);
   CHECK(!uno_optimize(NULL, model));
   CHECK(!uno_optimize(solver, NULL));
   // neither objective nor constraints
   CHECK(!uno_optimize(solver, model));
   uno_destroy_model(model);
   uno_destroy_solver(solver);
}

int main(void) {
   test_setters();
   test_optimize_with_invalid_arguments();
   if (number_failures > 0) {
      fprintf(stderr, "%d check(s) failed\n", number_failures);
      return EXIT_FAILURE;
   }
   printf("All option tests passed.\n");
   return EXIT_SUCCESS;
}
