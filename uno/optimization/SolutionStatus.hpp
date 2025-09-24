// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SOLUTIONSTATUS_H
#define UNO_SOLUTIONSTATUS_H

#include <string>

namespace uno {
   enum class SolutionStatus {
      NOT_OPTIMAL = 0,
      FEASIBLE_KKT_POINT, /* feasible stationary point */
      FEASIBLE_FJ_POINT, /* stationary point without constraint qualification */
      INFEASIBLE_STATIONARY_POINT, /* infeasible stationary point of constraint violation */
      FEASIBLE_SMALL_STEP,
      INFEASIBLE_SMALL_STEP,
      UNBOUNDED
   };

   std::string solution_status_to_message(SolutionStatus status);
} // namespace

#endif // UNO_SOLUTIONSTATUS_H