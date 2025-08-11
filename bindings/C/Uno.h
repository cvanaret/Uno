// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_C_API
#define UNO_C_API

#ifdef __cplusplus
extern "C" {
#endif

   // objective_function
   // - takes as input a vector "x" of size "number_variables" and an object "user_data", and
   // stores the objective value of "x" in "objective_value"
   // - returns an integer that is 0 if the evaluation succeeded, and positive otherwise
   int objective_function(const double* x, int number_variables, double* objective_value, void* user_data);

#ifdef __cplusplus
}
#endif

#endif // UNO_C_API