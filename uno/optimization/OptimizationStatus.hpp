// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIMIZATIONSTATUS_H
#define UNO_OPTIMIZATIONSTATUS_H

#include <string>

namespace uno {
   enum class OptimizationStatus {
      SUCCESS = 0,
      ITERATION_LIMIT,
      TIME_LIMIT,
      EVALUATION_ERROR,
      ALGORITHMIC_ERROR
   };

   std::string optimization_status_to_message(OptimizationStatus status);
} // namespace

#endif // UNO_OPTIMIZATIONSTATUS_H
