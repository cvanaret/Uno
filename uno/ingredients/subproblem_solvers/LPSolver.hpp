// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVER_H
#define UNO_LPSOLVER_H

#include "SubproblemSolver.hpp"

namespace uno {
   // forward declaration
   class Multipliers;

   class LPSolver: public SubproblemSolver {
   public:
      LPSolver() = default;
      ~LPSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override = 0;

      void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius, const Vector<double>& initial_point,
         Direction& direction, Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override = 0;

      [[nodiscard]] SolverWorkspace& get_workspace() override = 0;

   protected:
      static void compute_dual_displacements(const Subproblem& subproblem, Multipliers& direction_multipliers);
   };
} // namespace

#endif // UNO_LPSOLVER_H