// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAGRANGENEWTONSUBPROBLEM_H
#define UNO_LAGRANGENEWTONSUBPROBLEM_H

#include <memory>

namespace uno {
   // forward declarations
   class HessianModel;
   class Iterate;
   class Multipliers;
   class OptimizationProblem;
   class Options;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   class WarmstartInformation;

   // this class generalizes QP, LP (depends on the Hessian model) and EQP (as in interior-point methods)
   class LagrangeNewtonSubproblem {
   public:
      LagrangeNewtonSubproblem(const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
            bool use_regularization, double trust_region_radius, const Options& options);
      virtual ~LagrangeNewtonSubproblem();

      [[nodiscard]] const SymmetricMatrix<size_t, double>& get_lagrangian_hessian() const;

   protected:
      const OptimizationProblem& problem;
      Iterate& current_iterate;
      const Multipliers& current_multipliers;
      const std::unique_ptr<HessianModel> hessian_model;
      const double trust_region_radius;
   };
} // namespace

#endif // UNO_LAGRANGENEWTONSUBPROBLEM_H