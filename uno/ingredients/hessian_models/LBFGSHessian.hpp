// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HessianModel.hpp"

#ifndef UNO_LBFGSHESSIAN_H
#define UNO_LBFGSHESSIAN_H

namespace uno {
   class LBFGSHessian : public HessianModel {
   public:
      LBFGSHessian();
      ~LBFGSHessian() override = default;

      void initialize_statistics(Statistics& statistics, const Options& options) const override;
      void evaluate_hessian(Statistics& statistics, const OptimizationProblem& problem, const Vector<double>& primal_variables,
            const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
      void compute_hessian_vector_product(const OptimizationProblem& problem, const Vector<double>& vector,
         const Vector<double>& constraint_multipliers, Vector<double>& result) override;
   };
} // namespace

#endif // UNO_LBFGSHESSIAN_H