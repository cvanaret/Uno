// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EVALUATIONCACHE_H
#define UNO_EVALUATIONCACHE_H

#include "Evaluations.hpp"
#include "linear_algebra/COOSparsity.hpp"

namespace uno {
   // forward declaration
   class Model;

   class EvaluationCache {
   public:
      const size_t number_jacobian_nonzeros;
      COOSparsity jacobian_sparsity;
      // evaluations at current and trial iterates
      Evaluations current_evaluations;
      Evaluations trial_evaluations;

      explicit EvaluationCache(const Model& model);
   };
} // namespace

#endif // UNO_EVALUATIONCACHE_H