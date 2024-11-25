// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "OptimizationStatus.hpp"

namespace uno {
   std::string optimization_status_to_message(OptimizationStatus status) {
      if (status == OptimizationStatus::SUCCESS) {
         return "Success";
      }
      else if (status == OptimizationStatus::ITERATION_LIMIT) {
         return "Iteration limit";
      }
      else if (status == OptimizationStatus::TIME_LIMIT) {
         return "Time limit";
      }
      else if (status == OptimizationStatus::EVALUATION_ERROR) {
         return "Evaluation error";
      }
      else if (status == OptimizationStatus::ALGORITHMIC_ERROR) {
         return "Algorithmic error";
      }
      return "Unknown";
   }
} // namespace