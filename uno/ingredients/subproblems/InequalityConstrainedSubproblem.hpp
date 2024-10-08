// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYCONSTRAINEDSUBPROBLEM_H
#define UNO_INEQUALITYCONSTRAINEDSUBPROBLEM_H

#include "LagrangeNewtonSubproblem.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class Vector;

   class InequalityConstrainedSubproblem: LagrangeNewtonSubproblem {
   public:
      InequalityConstrainedSubproblem(const OptimizationProblem& problem, const Iterate& current_iterate, size_t number_variables,
            size_t number_hessian_nonzeros, bool use_regularization, double trust_region_radius, const Options& options):
            LagrangeNewtonSubproblem(problem, current_iterate, number_variables, number_hessian_nonzeros, use_regularization, trust_region_radius,
                  options) { }

      void variables_bounds(Vector<double>& variables_lower_bounds, Vector<double>& variables_upper_bounds) const;
      void linearized_constraint_bounds(const Vector<double>& constraints, Vector<double>& constraints_lower_bounds,
            Vector<double>& constraints_upper_bounds) const;

   protected:

   };
} // namespace

#endif // UNO_INEQUALITYCONSTRAINEDSUBPROBLEM_H