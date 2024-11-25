// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ITERATESTATUS_H
#define UNO_ITERATESTATUS_H

#include <string>

namespace uno {
   enum class IterateStatus {
      NOT_OPTIMAL = 0,
      FEASIBLE_KKT_POINT, /* feasible stationary point */
      FEASIBLE_FJ_POINT, /* stationary point without constraint qualification */
      INFEASIBLE_STATIONARY_POINT, /* infeasible stationary point of constraint violation */
      FEASIBLE_SMALL_STEP,
      INFEASIBLE_SMALL_STEP,
      UNBOUNDED
   };

   std::string iterate_status_to_message(IterateStatus status);
} // namespace

#endif // UNO_ITERATESTATUS_H
