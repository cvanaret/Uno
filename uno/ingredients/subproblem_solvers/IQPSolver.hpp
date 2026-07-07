// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_IQPSOLVER_H
#define UNO_IQPSOLVER_H

#include <memory>
#include "SubproblemSolver.hpp"

namespace uno {
   // forward declarations
   class LPSolver;

   // Inequality-constrained QP subproblem solver.
   //
   // IQPSolver is the inequality analogue of EQPSolver: it is the SubproblemSolver that takes a
   // Subproblem, builds a QuadraticProgram from it, and delegates the numerical solve to a backend
   // (BQPD/HiGHS) that knows only the QuadraticProgram. This mirrors EQPSolver delegating to a
   // DirectSymmetricIndefiniteLinearSolver that knows only a LinearSystem. The backend is held as an
   // LPSolver so that both curvature-free (LP) and quadratic (QP) subproblems can be wrapped; the
   // factory hands IQPSolver an LPSolver or a QPSolver as appropriate.
   class IQPSolver : public SubproblemSolver {
   public:
      explicit IQPSolver(std::unique_ptr<LPSolver> qp_solver);
      ~IQPSolver() override;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& initial_point, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] bool has_second_order_corrections() const override;
      void compute_second_order_correction(const Subproblem& subproblem, Direction& direction, Evaluations& trial_evaluations) override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   private:
      std::unique_ptr<LPSolver> qp_solver;

      static void compute_dual_direction(const Subproblem& subproblem, Multipliers& direction_multipliers);
   };
} // namespace

#endif // UNO_IQPSOLVER_H
