// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TRUSTREGIONSOLVER_H
#define UNO_TRUSTREGIONSOLVER_H

#include "InequalityConstrainedSolver.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class Statistics;
   class Subproblem;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class BoundConstrainedSolver: public InequalityConstrainedSolver {
   public:
      BoundConstrainedSolver() = default;
      ~BoundConstrainedSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override = 0;

      void solve(Statistics& statistics, Subproblem& subproblem, double trust_region_radius, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) override = 0;

      [[nodiscard]] virtual EvaluationSpace& get_evaluation_space() = 0;
   };
} // namespace

#endif // UNO_TRUSTREGIONSOLVER_H