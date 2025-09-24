// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BARRIERPROBLEM_H
#define UNO_BARRIERPROBLEM_H

#include "optimization/OptimizationProblem.hpp"

namespace uno {
   class BarrierProblem : public OptimizationProblem {
   public:
      explicit BarrierProblem(const Model& model, size_t number_variables, size_t number_constraints);
      ~BarrierProblem() override = default;

      virtual void generate_initial_iterate(Iterate& initial_iterate) const = 0;
      [[nodiscard]] virtual double compute_centrality_error(const Vector<double>& primals, const Multipliers& multipliers,
         double shift) const = 0;
   };

   inline BarrierProblem::BarrierProblem(const Model& model, size_t number_variables, size_t number_constraints):
      OptimizationProblem(model, number_variables, number_constraints) {
   }
} // namespace

#endif // UNO_BARRIERPROBLEM_H