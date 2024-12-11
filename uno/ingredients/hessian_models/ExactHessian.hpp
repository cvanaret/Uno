// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HessianModel.hpp"

namespace uno {
   // exact Hessian
   class ExactHessian : public HessianModel {
   public:
      ExactHessian();

      void initialize_statistics(Statistics& statistics, const Options& options) const override;
      void evaluate(Statistics& statistics, const OptimizationProblem& problem, const Vector<double>& primal_variables,
            const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
   };
} // namespace