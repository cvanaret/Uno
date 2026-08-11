// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEMSOLVER_H
#define UNO_SUBPROBLEMSOLVER_H

namespace uno {
   // forward declarations
   class Direction;
   class Evaluations;
   class Iterate;
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
      virtual void generate_initial_iterate(const Subproblem& subproblem, Iterate& initial_iterate, Evaluations& evaluations) = 0;
      virtual void compute_least_squares_multipliers(const Subproblem& subproblem, Iterate& iterate, Evaluations& evaluations) = 0;

      [[nodiscard]] virtual Direction& solve(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& initial_point, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) = 0;

      [[nodiscard]] virtual bool has_second_order_corrections() const = 0;
      virtual const Direction& compute_second_order_correction(const Subproblem& subproblem, const Iterate& current_iterate,
         const Vector<double>& constraints_SOC) = 0;

      [[nodiscard]] virtual const SolverWorkspace& get_workspace() const = 0;
   };
} // namespace

#endif // UNO_SUBPROBLEMSOLVER_H