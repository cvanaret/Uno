// Copyright (c) 2024-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "SolutionStatus.hpp"

namespace uno {
   std::string solution_status_to_message(SolutionStatus status) {
      if (status == SolutionStatus::FEASIBLE_KKT_POINT) {
         return "Feasible KKT point";
      }
      else if (status == SolutionStatus::FEASIBLE_FJ_POINT) {
         return "Feasible FJ point";
      }
      else if (status == SolutionStatus::INFEASIBLE_STATIONARY_POINT) {
         return "Infeasible stationary point";
      }
      else if (status == SolutionStatus::FEASIBLE_SMALL_STEP) {
         return "Feasible small step";
      }
      else if (status == SolutionStatus::INFEASIBLE_SMALL_STEP) {
         return "Infeasible small step";
      }
      else if (status == SolutionStatus::UNBOUNDED) {
         return "Unbounded problem";
      }
      return "Suboptimal point";
   }
} // namespace