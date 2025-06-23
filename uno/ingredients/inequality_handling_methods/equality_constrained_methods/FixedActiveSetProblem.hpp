// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FIXEDACTIVESETPROBLEM_H
#define UNO_FIXEDACTIVESETPROBLEM_H

#include "optimization/OptimizationProblem.hpp"

namespace uno {
   // forward declaration
   struct ActiveSet;

   class FixedActiveSetProblem: public OptimizationProblem {
   public:
      FixedActiveSetProblem(const OptimizationProblem& problem, const Vector<double>& current_primals, const ActiveSet& active_set,
         double trust_region_radius);

      [[nodiscard]] virtual double variable_lower_bound(size_t variable_index) const;
      [[nodiscard]] virtual double variable_upper_bound(size_t variable_index) const;
      [[nodiscard]] virtual double constraint_lower_bound(size_t constraint_index) const;
      [[nodiscard]] virtual double constraint_upper_bound(size_t constraint_index) const;

   protected:
      const OptimizationProblem& problem;
      const Vector<double>& direction_primals;
      const ActiveSet& active_set;
      double trust_region_radius;

      const double activity_tolerance{1e-6}; // TODO option
      Vector<double> variables_lower_bounds;
      Vector<double> variables_upper_bounds;
      Vector<double> constraints_lower_bounds;
      Vector<double> constraints_upper_bounds;
   };
} // namespace

#endif // UNO_FIXEDACTIVESETPROBLEM_H