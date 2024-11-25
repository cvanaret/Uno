// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "IterateStatus.hpp"

namespace uno {
   std::string iterate_status_to_message(IterateStatus status) {
      if (status == IterateStatus::FEASIBLE_KKT_POINT) {
         return "Feasible KKT point";
      }
      else if (status == IterateStatus::FEASIBLE_FJ_POINT) {
         return "Feasible FJ point";
      }
      else if (status == IterateStatus::INFEASIBLE_STATIONARY_POINT) {
         return "Infeasible stationary point";
      }
      else if (status == IterateStatus::FEASIBLE_SMALL_STEP) {
         return "Feasible small step";
      }
      else if (status == IterateStatus::INFEASIBLE_SMALL_STEP) {
         return "Infeasible small step";
      }
      else if (status == IterateStatus::UNBOUNDED) {
         return "Unbounded problem";
      }
      return "Suboptimal point";
   }
} // namespace