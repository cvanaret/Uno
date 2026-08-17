// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HessianBuffer.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "tools/Logger.hpp"

namespace uno {
   void HessianBuffer::initialize(size_t number_hessian_nonzeros) {
      this->hessian_buffer.resize(number_hessian_nonzeros);
   }

   void HessianBuffer::evaluate_hessian(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         View<double> hessian_values) {
      assert(this->hessian_buffer.size() == hessian_values.size());

      // try to evaluate the Hessian. Upon failure, keep the previous one
      try {
         subproblem.evaluate_lagrangian_hessian(statistics, current_iterate, this->hessian_buffer.view());
         // success: copy the new Hessian into the linear system
         hessian_values = this->hessian_buffer;
      }
      catch (const HessianEvaluationError&) {
         if (this->first_evaluation) {
            hessian_values.fill(1.);
            DEBUG << "The initial Hessian could not be evaluated, using the identity\n";
         }
         else {
            // keep the existing Hessian
            DEBUG << "The Hessian could not be evaluated, using the previous one\n";
         }
      }
      this->first_evaluation = false;
   }
} // namespace