// Copyright (c) 2025 Charlie Vanaret, Alexis Montoison
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTANTS_H
#define UNO_CONSTANTS_H

// Optimization sense
const int32_t UNO_MINIMIZE = 1;
const int32_t UNO_MAXIMIZE = -1;

// Lagrange multiplier sign convention
const double UNO_MULTIPLIER_POSITIVE =  1.0;
const double UNO_MULTIPLIER_NEGATIVE = -1.0;

// Problem type: 'L' = Linear, 'Q' = Quadratic, 'N' = Nonlinear
const char UNO_PROBLEM_LINEAR    = 'L';
const char UNO_PROBLEM_QUADRATIC = 'Q';
const char UNO_PROBLEM_NONLINEAR = 'N';

// Base indexing style: 0-based (C) or 1-based (Fortran)
const int32_t UNO_ZERO_BASED_INDEXING = 0;
const int32_t UNO_ONE_BASED_INDEXING  = 1;

// Triangular part: 'L' = lower, 'U' = upper
const char UNO_LOWER_TRIANGLE = 'L';
const char UNO_UPPER_TRIANGLE = 'U';

// Optimization status
const int32_t UNO_SUCCESS = 0;
const int32_t UNO_ITERATION_LIMIT = 1;
const int32_t UNO_TIME_LIMIT = 2;
const int32_t UNO_EVALUATION_ERROR = 3;
const int32_t UNO_ALGORITHMIC_ERROR = 4;

// Iterate status
const int32_t UNO_NOT_OPTIMAL = 0;
const int32_t UNO_FEASIBLE_KKT_POINT = 1; // feasible stationary point
const int32_t UNO_FEASIBLE_FJ_POINT = 2; // stationary point without constraint qualification
const int32_t UNO_INFEASIBLE_STATIONARY_POINT = 3; // infeasible stationary point of constraint violation
const int32_t UNO_FEASIBLE_SMALL_STEP = 4;
const int32_t UNO_INFEASIBLE_SMALL_STEP = 5;
const int32_t UNO_UNBOUNDED = 6;

// current Uno version is 2.0.3
const int32_t UNO_VERSION_MAJOR = 2;
const int32_t UNO_VERSION_MINOR = 0;
const int32_t UNO_VERSION_PATCH = 3;

#endif // UNO_CONSTANTS_H