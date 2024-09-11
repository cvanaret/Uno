// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEMSTATUS_H
#define UNO_SUBPROBLEMSTATUS_H

namespace uno {
   enum class SubproblemStatus {
      OPTIMAL = 0,
      UNBOUNDED_PROBLEM,
      INFEASIBLE,
      ERROR
   };
} // namespace

#endif // UNO_SUBPROBLEMSTATUS_H
