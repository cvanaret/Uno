// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EQPSOLVER_H
#define UNO_EQPSOLVER_H

#include <memory>
#include "SubproblemSolver.hpp"
#include "optimization/Direction.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class Options;

   class EQPSolver: public SubproblemSolver {
   public:
      explicit EQPSolver(const Options& options);
      ~EQPSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;
      void generate_initial_iterate(const Subproblem& subproblem, Iterate& initial_iterate, Evaluations& evaluations) override;
      void compute_least_squares_multipliers(const Subproblem& subproblem, Iterate& iterate, Evaluations& evaluations,
         double multipliers_threshold) override;

      [[nodiscard]] const Direction& solve(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& initial_point, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] bool has_second_order_corrections() const override;
      void initialize_second_order_corrections(const Subproblem& subproblem, const Iterate& current_iterate,
         const Iterate& trial_iterate, Evaluations& current_evaluations, Evaluations& trial_evaluations) override;
      const Direction& compute_second_order_correction(const Subproblem& subproblem, const Iterate& current_iterate) override;
      void update_second_order_corrections(const Subproblem& subproblem, const Iterate& trial_iterate,
         Evaluations& trial_evaluations) override;

      [[nodiscard]] const SolverWorkspace& get_workspace() const override;

   protected:
      Direction direction;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
      bool analysis_performed{false};

      bool SOC_initialized{false};
      Vector<double> constraints_SOC;
      Vector<double> constraints_buffer_SOC;
      Direction direction_SOC;
   };
} // namespace

#endif // UNO_EQPSOLVER_H