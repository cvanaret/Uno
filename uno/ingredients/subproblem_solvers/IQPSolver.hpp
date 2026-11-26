// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_IQPSOLVER_H
#define UNO_IQPSOLVER_H

#include <memory>
#include "SubproblemSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"

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
   // forward declaration
   class Evaluations;

   class IQPSolver : public SubproblemSolver {
   public:
      explicit IQPSolver(std::unique_ptr<LPSolver> qp_solver);
      ~IQPSolver() override;

      void initialize_memory(const Subproblem& subproblem) override;
      [[nodiscard]] bool requires_hessian_matrix(const Subproblem& subproblem) const override;
      void compute_least_squares_multipliers(const Subproblem& subproblem, Iterate& iterate, Evaluations& evaluations,
         double multipliers_threshold) override;

      [[nodiscard]] const Direction& solve(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& initial_point, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] bool has_second_order_corrections() const override;
      void initialize_second_order_corrections(const Subproblem& subproblem, const Iterate& current_iterate,
         const Iterate& trial_iterate, Evaluations& current_evaluations, Evaluations& trial_evaluations) override;
      const Direction& compute_second_order_correction(Statistics& statistics, const Subproblem& subproblem,
         const Iterate& current_iterate) override;
      void update_second_order_corrections(const Subproblem& subproblem, const Iterate& trial_iterate,
         Evaluations& trial_evaluations) override;

      [[nodiscard]] const SolverWorkspace& get_workspace() const override;

   private:
      Direction direction;
      std::unique_ptr<LPSolver> qp_solver;
      // second-order corrections: constants c(x_k + d) - J_k d of the linearized constraints
      Direction direction_SOC;
      Vector<double> constraints_SOC{};
      Vector<double> jacobian_vector_product_SOC{};
      Vector<double> initial_point_SOC{};
      Evaluations* current_evaluations_SOC{nullptr};

      void compute_SOC_constraint_constants(const Subproblem& subproblem, const Iterate& trial_iterate,
         const Direction& step, Evaluations& trial_evaluations);

      static void compute_dual_direction(const Subproblem& subproblem, const Iterate& current_iterate, Multipliers& direction_multipliers);
   };
} // namespace

#endif // UNO_IQPSOLVER_H
