// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVER_H
#define UNO_LPSOLVER_H

#include "InequalityConstrainedSolver.hpp"

namespace uno {
   class LPSolver: public InequalityConstrainedSolver {
   public:
      LPSolver() = default;
      ~LPSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override = 0;

      void solve(Statistics& statistics, Subproblem& subproblem, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) override = 0;

      [[nodiscard]] EvaluationSpace& get_evaluation_space() override = 0;
   };
} // namespace

#endif // UNO_LPSOLVER_H