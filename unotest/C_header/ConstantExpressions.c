// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// The integer constants of the C API are constant expressions in C: they can be used as case labels, and the version
// numbers can be tested by the preprocessor (compile-time checks)

#include "../../interfaces/C/Uno_C_API.h"

#if !defined(UNO_VERSION_MAJOR) || !defined(UNO_VERSION_MINOR) || !defined(UNO_VERSION_PATCH)
#error "the version numbers must be preprocessor macros"
#endif
#if UNO_VERSION_MAJOR < 0 || UNO_VERSION_MINOR < 0 || UNO_VERSION_PATCH < 0
#error "invalid version numbers"
#endif

const char* uno_header_optimization_status_name(uno_int status) {
   switch (status) {
      case UNO_SUCCESS: return "success";
      case UNO_ITERATION_LIMIT: return "iteration limit";
      case UNO_TIME_LIMIT: return "time limit";
      case UNO_EVALUATION_ERROR: return "evaluation error";
      case UNO_ALGORITHMIC_ERROR: return "algorithmic error";
      case UNO_USER_TERMINATION: return "user termination";
      case UNO_USER_INTERRUPTION: return "user interruption";
      default: return "unknown";
   }
}

const char* uno_header_solution_status_name(uno_int status) {
   switch (status) {
      case UNO_NOT_OPTIMAL: return "not optimal";
      case UNO_FEASIBLE_KKT_POINT: return "feasible KKT point";
      case UNO_FEASIBLE_FJ_POINT: return "feasible FJ point";
      case UNO_INFEASIBLE_STATIONARY_POINT: return "infeasible stationary point";
      case UNO_FEASIBLE_SMALL_STEP: return "feasible small step";
      case UNO_INFEASIBLE_SMALL_STEP: return "infeasible small step";
      case UNO_DIVERGING_ITERATE: return "diverging iterate";
      case UNO_UNBOUNDED_OBJECTIVE: return "unbounded objective";
      default: return "unknown";
   }
}

int uno_header_is_success(uno_int status) {
   switch (status) {
      case UNO_OPTIONS_SUCCESS: return 1; // all the success codes are 0
      case UNO_OPTIONS_INVALID_ARGUMENT:
      case UNO_OPTIONS_NOT_FOUND:
      case UNO_OPTIONS_WRONG_TYPE:
      case UNO_OPTIONS_INVALID_VALUE:
      case UNO_OPTIONS_INVALID_FILE:
      case UNO_OPTIONS_ERROR:
      case UNO_OPTIONS_NOT_SET:
         return 0;
      default: return -1;
   }
}

int uno_header_other_constants(uno_int value, char triangle) {
   int sum = 0;
   switch (value) {
      case UNO_MINIMIZE: sum += 1; break; // UNO_MULTIPLIER_POSITIVE has the same value
      case UNO_MAXIMIZE: sum += 2; break;
      default: break;
   }
   switch (value) {
      case UNO_ZERO_BASED_INDEXING: sum += 4; break;
      case UNO_ONE_BASED_INDEXING: sum += 8; break;
      default: break;
   }
   switch (value) {
      case UNO_OPTION_TYPE_INTEGER: case UNO_OPTION_TYPE_DOUBLE: case UNO_OPTION_TYPE_BOOL: case UNO_OPTION_TYPE_STRING:
      case UNO_OPTION_TYPE_NOT_FOUND: sum += 16; break;
      default: break;
   }
   switch (value) {
      case UNO_RESULT_SUCCESS: case UNO_RESULT_INVALID_ARGUMENT: case UNO_RESULT_NOT_AVAILABLE:
      case UNO_RESULT_INVALID_INDEX: case UNO_RESULT_ERROR: sum += 32; break;
      default: break;
   }
   switch (value) {
      case UNO_OPTIMIZE_SUCCESS: case UNO_OPTIMIZE_INVALID_ARGUMENT: case UNO_OPTIMIZE_INVALID_MODEL:
      case UNO_OPTIMIZE_INVALID_OPTIONS: case UNO_OPTIMIZE_ERROR: sum += 64; break;
      default: break;
   }
   switch (value) {
      case UNO_MODEL_SUCCESS: case UNO_MODEL_INVALID_ARGUMENT: case UNO_MODEL_INVALID_INDEX:
      case UNO_MODEL_INVALID_VALUE: case UNO_MODEL_ERROR: sum += 128; break;
      default: break;
   }
   switch (value) {
      case UNO_SOLVER_SUCCESS: case UNO_SOLVER_INVALID_ARGUMENT: case UNO_SOLVER_ERROR: sum += 256; break;
      default: break;
   }
   switch (value) {
      case UNO_MULTIPLIER_NEGATIVE: sum += 512; break;
      default: break;
   }
   switch (triangle) {
      case UNO_LOWER_TRIANGLE: sum += 1024; break;
      case UNO_UPPER_TRIANGLE: sum += 2048; break;
      default: break;
   }
   return sum;
}
