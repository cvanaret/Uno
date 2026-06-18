// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVER_H
#define UNO_LPSOLVER_H

namespace uno {
   // forward declarations
   class Direction;
   class QuadraticProgram;
   class SolverWorkspace;
   class Statistics;
   class Subproblem;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   // Numerical backend that solves a (Subproblem-free) QuadraticProgram. An LPSolver is the base
   // capability: it is guaranteed to solve curvature-free programs (LPs). BQPD and HiGHS implement
   // the QPSolver refinement (QPSolver : LPSolver), so they are valid LPSolvers too; a simplex-only
   // backend would implement LPSolver directly.
   //
   // The backend owns its native QuadraticProgram (only the backend knows the format). IQPSolver
   // retrieves it via get_quadratic_program() and populates it from a Subproblem before each solve;
   // solve() itself never sees a Subproblem.
   class LPSolver {
   public:
      LPSolver() = default;
      virtual ~LPSolver();

      // allocate native storage; also constructs and initializes the owned QuadraticProgram
      virtual void initialize_memory(const Subproblem& subproblem) = 0;

      // the QuadraticProgram owned by this backend, populated by IQPSolver before each solve
      [[nodiscard]] virtual QuadraticProgram& get_quadratic_program() = 0;

      // solve the (already built) QuadraticProgram. Fills direction.primals, direction.multipliers
      // (raw subproblem multipliers, before dual-displacement mapping), direction.status and
      // direction.subproblem_objective.
      virtual void solve(Statistics& statistics, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) = 0;

      [[nodiscard]] virtual SolverWorkspace& get_workspace() = 0;
   };
} // namespace

#endif // UNO_LPSOLVER_H
