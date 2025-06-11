// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EQPPROBLEM_H
#define UNO_EQPPROBLEM_H

#include "optimization/OptimizationProblem.hpp"

namespace uno {
   // forward declaration
   class ActiveSet;

   class EQPProblem: public OptimizationProblem {
   public:
      EQPProblem(const OptimizationProblem& problem, const ActiveSet& active_set);

   protected:
      const OptimizationProblem& first_reformulation;
      const ActiveSet& active_set;
   };
} // namespace

#endif // UNO_EQPPROBLEM_H