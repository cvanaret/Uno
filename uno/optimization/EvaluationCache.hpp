// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EVALUATIONCACHE_H
#define UNO_EVALUATIONCACHE_H

#include "Evaluations.hpp"

namespace uno {
   class EvaluationCache {
   public:
      // evaluations at current and trial iterates
      Evaluations current_evaluations;
      Evaluations trial_evaluations;

      EvaluationCache(size_t number_variables, size_t number_constraints);
   };
} // namespace

#endif // UNO_EVALUATIONCACHE_H