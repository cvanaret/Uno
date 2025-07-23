// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSOLVER_H
#define UNO_QPSOLVER_H

#include "LPSolver.hpp"

namespace uno {
   class QPSolver : public LPSolver {
   public:
      QPSolver() = default;
      ~QPSolver() override = default;

      void initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model,
         const RegularizationStrategy<double>& regularization_strategy) override = 0;

      void solve_inequality_constrained_subproblem(Statistics& statistics, Subproblem& subproblem,
         const Vector<double>& initial_point, Direction& direction, const WarmstartInformation& warmstart_information) override = 0;

      [[nodiscard]] double hessian_quadratic_product(const Vector<double>& vector) const override = 0;
   };
} // namespace

#endif // UNO_QPSOLVER_H