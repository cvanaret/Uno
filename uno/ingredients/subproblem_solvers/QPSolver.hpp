// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSOLVER_H
#define UNO_QPSOLVER_H

#include <vector>
#include "LPSolver.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class Iterate;
   class HessianModel;
   class OptimizationProblem;
   class Options;
   template <typename ElementType>
   class RectangularMatrix;
   template <typename ElementType>
   class SparseVector;
   class Statistics;
   struct WarmstartInformation;

   class QPSolver : public LPSolver {
   public:
      QPSolver(): LPSolver() { }
      ~QPSolver() override = default;

      void solve_LP(const OptimizationProblem& problem, Iterate& current_iterate, const Vector<double>& initial_point, Direction& direction,
            double trust_region_radius, const WarmstartInformation& warmstart_information) override = 0;

      virtual void solve_QP(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate, const Vector<double>& current_multipliers,
            const Vector<double>& initial_point, Direction& direction, HessianModel& hessian_model, double trust_region_radius,
            const WarmstartInformation& warmstart_information) = 0;
   };
} // namespace

#endif // UNO_QPSOLVER_H