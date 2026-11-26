// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "../interfaces/C/Uno_C_API.h"

namespace py = pybind11;

namespace uno {
   // the integer constants of the C API are enumerators: they are cast, since pybind11 cannot convert an unregistered
   // enum type
   void define_Constants(py::module& module) {
      // Optimization sense
      module.attr("MINIMIZE") = static_cast<uno_int>(UNO_MINIMIZE);
      module.attr("MAXIMIZE") = static_cast<uno_int>(UNO_MAXIMIZE);

      // Lagrange multiplier sign convention
      module.attr("MULTIPLIER_POSITIVE") = static_cast<uno_int>(UNO_MULTIPLIER_POSITIVE);
      module.attr("MULTIPLIER_NEGATIVE") = static_cast<uno_int>(UNO_MULTIPLIER_NEGATIVE);

      // Problem type: "LP" = linear, "QP" = quadratic, "NLP" = nonlinear
      module.attr("PROBLEM_LINEAR") = UNO_PROBLEM_LINEAR;
      module.attr("PROBLEM_QUADRATIC") = UNO_PROBLEM_QUADRATIC;
      module.attr("PROBLEM_NONLINEAR") = UNO_PROBLEM_NONLINEAR;

      // Base indexing style: 0-based (C) or 1-based (Fortran)
      module.attr("ZERO_BASED_INDEXING") = static_cast<uno_int>(UNO_ZERO_BASED_INDEXING);
      module.attr("ONE_BASED_INDEXING") = static_cast<uno_int>(UNO_ONE_BASED_INDEXING);

      // Triangular part: 'L' = lower, 'U' = upper
      module.attr("LOWER_TRIANGLE") = static_cast<char>(UNO_LOWER_TRIANGLE);
      module.attr("UPPER_TRIANGLE") = static_cast<char>(UNO_UPPER_TRIANGLE);

      // Optimization status
      module.attr("SUCCESS") = static_cast<uno_int>(UNO_SUCCESS);
      module.attr("ITERATION_LIMIT") = static_cast<uno_int>(UNO_ITERATION_LIMIT);
      module.attr("TIME_LIMIT") = static_cast<uno_int>(UNO_TIME_LIMIT);
      module.attr("EVALUATION_ERROR") = static_cast<uno_int>(UNO_EVALUATION_ERROR);
      module.attr("ALGORITHMIC_ERROR") = static_cast<uno_int>(UNO_ALGORITHMIC_ERROR);
      module.attr("USER_TERMINATION") = static_cast<uno_int>(UNO_USER_TERMINATION);
      module.attr("USER_INTERRUPTION") = static_cast<uno_int>(UNO_USER_INTERRUPTION);

      // Iterate status
      module.attr("NOT_OPTIMAL") = static_cast<uno_int>(UNO_NOT_OPTIMAL);
      module.attr("FEASIBLE_KKT_POINT") = static_cast<uno_int>(UNO_FEASIBLE_KKT_POINT);
      module.attr("FEASIBLE_FJ_POINT") = static_cast<uno_int>(UNO_FEASIBLE_FJ_POINT);
      module.attr("INFEASIBLE_STATIONARY_POINT") = static_cast<uno_int>(UNO_INFEASIBLE_STATIONARY_POINT);
      module.attr("FEASIBLE_SMALL_STEP") = static_cast<uno_int>(UNO_FEASIBLE_SMALL_STEP);
      module.attr("INFEASIBLE_SMALL_STEP") = static_cast<uno_int>(UNO_INFEASIBLE_SMALL_STEP);
      module.attr("DIVERGING_ITERATE") = static_cast<uno_int>(UNO_DIVERGING_ITERATE);
      module.attr("UNBOUNDED_OBJECTIVE") = static_cast<uno_int>(UNO_UNBOUNDED_OBJECTIVE);

      // current Uno version
      module.attr("VERSION_MAJOR") = static_cast<uno_int>(UNO_VERSION_MAJOR);
      module.attr("VERSION_MINOR") = static_cast<uno_int>(UNO_VERSION_MINOR);
      module.attr("VERSION_PATCH") = static_cast<uno_int>(UNO_VERSION_PATCH);
   }
} // namespace