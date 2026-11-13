// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "EvaluationCache.hpp"
#include "model/Model.hpp"

namespace uno {
   EvaluationCache::EvaluationCache(const Model& model):
         number_jacobian_nonzeros(model.number_jacobian_nonzeros()),
         current_evaluations(model),
         trial_evaluations(model) {
   }
} // namespace