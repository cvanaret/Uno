// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINTERIORPOINTSYSTEM_H
#define UNO_PRIMALDUALINTERIORPOINTSYSTEM_H

#include "LagrangeNewtonSubproblem.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class Vector;

   class PrimalDualInteriorPointSystem: LagrangeNewtonSubproblem {
   public:
      PrimalDualInteriorPointSystem(const OptimizationProblem& problem, const Iterate& current_iterate, size_t number_variables,
            size_t number_hessian_nonzeros, bool use_regularization, double trust_region_radius, const Options& options):
            LagrangeNewtonSubproblem(problem, current_iterate, number_variables, number_hessian_nonzeros, use_regularization, trust_region_radius,
                  options) { }

      void right_hand_side(Vector<double>& rhs) const;

   protected:

   };
} // namespace

#endif // UNO_PRIMALDUALINTERIORPOINTSYSTEM_H