// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NUMBERMODELEVALUATIONS_H
#define UNO_NUMBERMODELEVALUATIONS_H

#include <cstddef>

namespace uno {
   class NumberModelEvaluations {
   public:
      size_t objective{0};
      size_t constraints{0};
      size_t objective_gradient{0};
      size_t jacobian{0};
      size_t hessian{0};

      void reset() {
         this->objective = this->constraints = this->objective_gradient = this->jacobian = this->hessian = 0;
      }
   };
} // namespace

#endif //UNO_NUMBERMODELEVALUATIONS_H