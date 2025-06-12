// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FIXEDACTIVESETPROBLEM_H
#define UNO_FIXEDACTIVESETPROBLEM_H

#include "optimization/OptimizationProblem.hpp"

namespace uno {
   // forward declaration
   class ActiveSet;

   class FixedActiveSetProblem: public OptimizationProblem {
   public:
      FixedActiveSetProblem(const OptimizationProblem& problem, const ActiveSet& active_set);

   protected:
      const OptimizationProblem& first_reformulation;
      const ActiveSet& active_set;
   };
} // namespace

#endif // UNO_FIXEDACTIVESETPROBLEM_H