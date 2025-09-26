// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "../bindings/C/Uno_C_API.h"

namespace py = pybind11;

namespace uno {
   void define_Constants(py::module& module) {
      // Optimization sense
      module.attr("MINIMIZE") = UNO_MINIMIZE;
      module.attr("MAXIMIZE") = UNO_MAXIMIZE;

      // Lagrange multiplier sign convention
      module.attr("MULTIPLIER_POSITIVE") = UNO_MULTIPLIER_POSITIVE;
      module.attr("MULTIPLIER_NEGATIVE") = UNO_MULTIPLIER_NEGATIVE;

      // Problem type: 'L' = Linear, 'Q' = Quadratic, 'N' = Nonlinear
      module.attr("PROBLEM_LINEAR") = UNO_PROBLEM_LINEAR;
      module.attr("PROBLEM_QUADRATIC") = UNO_PROBLEM_QUADRATIC;
      module.attr("PROBLEM_NONLINEAR") = UNO_PROBLEM_NONLINEAR;

      // Base indexing style: 0-based (C) or 1-based (Fortran)
      module.attr("ZERO_BASED_INDEXING") = UNO_ZERO_BASED_INDEXING;
      module.attr("ONE_BASED_INDEXING") = UNO_ONE_BASED_INDEXING;

      // Triangular part: 'L' = lower, 'U' = upper
      module.attr("LOWER_TRIANGLE") = UNO_LOWER_TRIANGLE;
      module.attr("UPPER_TRIANGLE") = UNO_UPPER_TRIANGLE;

      // Optimization status
      module.attr("SUCCESS") = UNO_SUCCESS;
      module.attr("ITERATION_LIMIT") = UNO_ITERATION_LIMIT;
      module.attr("TIME_LIMIT") = UNO_TIME_LIMIT;
      module.attr("EVALUATION_ERROR") = UNO_EVALUATION_ERROR;
      module.attr("ALGORITHMIC_ERROR") = UNO_ALGORITHMIC_ERROR;

      // Iterate status
      module.attr("NOT_OPTIMAL") = UNO_NOT_OPTIMAL;
      module.attr("FEASIBLE_KKT_POINT") = UNO_FEASIBLE_KKT_POINT;
      module.attr("FEASIBLE_FJ_POINT") = UNO_FEASIBLE_FJ_POINT;
      module.attr("INFEASIBLE_STATIONARY_POINT") = UNO_INFEASIBLE_STATIONARY_POINT;
      module.attr("FEASIBLE_SMALL_STEP") = UNO_FEASIBLE_SMALL_STEP;
      module.attr("INFEASIBLE_SMALL_STEP") = UNO_INFEASIBLE_SMALL_STEP;
      module.attr("UNBOUNDED") = UNO_UNBOUNDED;

      // current Uno version is 2.0.3
      module.attr("VERSION_MAJOR") = UNO_VERSION_MAJOR;
      module.attr("VERSION_MINOR") = UNO_VERSION_MINOR;
      module.attr("VERSION_PATCH") = UNO_VERSION_PATCH;
   }
} // namespace