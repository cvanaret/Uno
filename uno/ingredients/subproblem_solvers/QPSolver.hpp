// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSOLVER_H
#define UNO_QPSOLVER_H

#include "LPSolver.hpp"

namespace uno {
   class QPSolver : public LPSolver {
   public:
      QPSolver(): LPSolver() { }
      ~QPSolver() override = default;

      void initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model,
         RegularizationStrategy<double>& regularization_strategy) override = 0;

      void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
         const Vector<double>& initial_point, Direction& direction, HessianModel& hessian_model,
         RegularizationStrategy<double>& regularization_strategy, double trust_region_radius,
         const WarmstartInformation& warmstart_information) override = 0;
   };
} // namespace

#endif // UNO_QPSOLVER_H