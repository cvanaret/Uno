#include "EvaluationCache.hpp"

namespace uno {
   EvaluationCache::EvaluationCache(size_t number_variables, size_t number_constraints):
      current_evaluations(number_variables, number_constraints),
      trial_evaluations(number_variables, number_constraints) {
   }
} // namespace