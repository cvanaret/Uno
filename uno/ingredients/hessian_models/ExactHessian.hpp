// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HessianModel.hpp"

namespace uno {
   // forward declaration
   class Options;

   // exact Hessian
   class ExactHessian : public HessianModel {
   public:
      ExactHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options);

      void evaluate(Statistics& statistics, const OptimizationProblem& problem, const Vector<double>& primal_variables,
            const Vector<double>& constraint_multipliers) override;
   };
} // namespace