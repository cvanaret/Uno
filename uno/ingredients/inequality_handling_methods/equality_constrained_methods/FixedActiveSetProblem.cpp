// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FixedActiveSetProblem.hpp"

namespace uno {
   FixedActiveSetProblem::FixedActiveSetProblem(const OptimizationProblem& problem, const ActiveSet &active_set):
         OptimizationProblem(problem.model), first_reformulation(problem), active_set(active_set) {
   }
} // namespace