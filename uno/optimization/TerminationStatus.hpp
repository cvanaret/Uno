// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TERMINATIONSTATUS_H
#define UNO_TERMINATIONSTATUS_H

enum TerminationStatus {
   NOT_OPTIMAL = 0,
   KKT_POINT, /* feasible stationary point */
   FJ_POINT, /* infeasible stationary point */
   FEASIBLE_SMALL_STEP,
   INFEASIBLE_SMALL_STEP
};

#endif // UNO_TERMINATIONSTATUS_H