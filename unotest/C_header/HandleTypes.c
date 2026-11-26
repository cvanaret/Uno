// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// The C API handles are typed in C too, and a handle converts implicitly from/to void* (C source compatibility)

#include "../../interfaces/C/Uno_C_API.h"

uno_int uno_header_handle_types(void) {
   uno_model* model = uno_create_unconstrained_model(UNO_PROBLEM_NONLINEAR, 1, UNO_ZERO_BASED_INDEXING);
   void* solver = uno_create_solver(); // implicit conversion uno_solver* -> void*
   const uno_solver* const_solver = (const uno_solver*)solver;
   uno_int status = -1;
   const uno_int result = uno_get_optimization_status(const_solver, &status);
   uno_destroy_solver(solver); // implicit conversion void* -> uno_solver*
   uno_destroy_model(model);
   return result;
}
