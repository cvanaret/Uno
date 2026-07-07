// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEMSOLVER_H
#define UNO_SUBPROBLEMSOLVER_H

namespace uno {
   // forward declarations
   class Direction;
   class Evaluations;
   class Multipliers;
   class Statistics;
   class SolverWorkspace;
   class Subproblem;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class SubproblemSolver {
   public:
      SubproblemSolver() = default;
      virtual ~SubproblemSolver() = default;

      virtual void initialize_memory(const Subproblem& subproblem) = 0;

      virtual void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& initial_point, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) = 0;

      [[nodiscard]] virtual bool has_second_order_corrections() const = 0;
      virtual void compute_second_order_correction(const Subproblem& subproblem, Direction& direction,
         Evaluations& trial_evaluations) = 0;

      [[nodiscard]] virtual SolverWorkspace& get_workspace() = 0;
   };
} // namespace

#endif // UNO_SUBPROBLEMSOLVER_H