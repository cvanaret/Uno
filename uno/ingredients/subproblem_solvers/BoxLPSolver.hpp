// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BOXLPSOLVER_H
#define UNO_BOXLPSOLVER_H

#include <vector>
#include "SubproblemSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "SolverWorkspace.hpp"

namespace uno {
   class BoxLPSolverWorkspace: public SolverWorkspace {
   public:
      BoxLPSolverWorkspace() = default;

      [[nodiscard]] double compute_hessian_quadratic_form(const Subproblem& /*subproblem*/, const Iterate& /*current_iterate*/,
            const Vector<double>& /*vector*/) const override {
         return 0.;
      }

      Vector<double> objective_gradient;
   };

   class BoxLPSolver: public SubproblemSolver {
   public:
      BoxLPSolver() = default;
      ~BoxLPSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;
      void generate_initial_iterate(const Subproblem& subproblem, Iterate& initial_iterate, Evaluations& evaluations) override;
      void compute_least_squares_multipliers(const Subproblem& subproblem, Iterate& iterate, Evaluations& evaluations,
         double multipliers_threshold) override;

      [[nodiscard]] Direction& solve(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& initial_point, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] bool has_second_order_corrections() const override;
      const Direction& compute_second_order_correction(const Subproblem& subproblem, const Iterate& current_iterate,
         const Vector<double>& constraints_SOC) override;

      [[nodiscard]] const SolverWorkspace& get_workspace() const override;

   protected:
      Direction direction;
      std::vector<double> variable_lower_bounds;
      std::vector<double> variable_upper_bounds;
      BoxLPSolverWorkspace workspace{};
   };
} // namespace

#endif // UNO_BOXLPSOLVER_H