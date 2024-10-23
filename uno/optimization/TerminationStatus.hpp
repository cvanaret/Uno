// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TERMINATIONSTATUS_H
#define UNO_TERMINATIONSTATUS_H

namespace uno {
   enum class TerminationStatus {
      NOT_OPTIMAL = 0,
      FEASIBLE_KKT_POINT, /* feasible stationary point */
      FEASIBLE_FJ_POINT, /* stationary point without constraint qualification */
      INFEASIBLE_STATIONARY_POINT, /* infeasible stationary point of constraint violation */
      FEASIBLE_SMALL_STEP,
      INFEASIBLE_SMALL_STEP,
      UNBOUNDED
   };

   inline std::string status_to_message(TerminationStatus status) {
      if (status == TerminationStatus::FEASIBLE_KKT_POINT) {
         return "Converged with feasible KKT point";
      }
      else if (status == TerminationStatus::FEASIBLE_FJ_POINT) {
         return "Converged with feasible FJ point";
      }
      else if (status == TerminationStatus::INFEASIBLE_STATIONARY_POINT) {
         return "Converged with infeasible stationary point";
      }
      else if (status == TerminationStatus::FEASIBLE_SMALL_STEP) {
         return "Terminated with feasible small step";
      }
      else if (status == TerminationStatus::INFEASIBLE_SMALL_STEP) {
         return "Failed with infeasible small step";
      }
      else if (status == TerminationStatus::UNBOUNDED) {
         return "Terminated with unbounded problem";
      }
      else {
         return "Failed with suboptimal point";
      }
   }
} // namespace

#endif // UNO_TERMINATIONSTATUS_H
