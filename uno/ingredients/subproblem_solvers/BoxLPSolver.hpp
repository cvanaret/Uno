// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BOXLPSOLVER_H
#define UNO_BOXLPSOLVER_H

#include <vector>
#include "InequalityConstrainedSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "SolverWorkspace.hpp"

namespace uno {
   class BoxLPSolverWorkspace: public SolverWorkspace {
   public:
      BoxLPSolverWorkspace() = default;

      void evaluate_jacobian(const OptimizationProblem& /*problem*/, const Vector<double>& /*primals*/) override { }
      [[nodiscard]] double compute_hessian_quadratic_product(const Subproblem& /*subproblem*/, const Vector<double>& /*vector*/) const override {
         return 0.;
      }

      Vector<double> objective_gradient;
   };

   class BoxLPSolver: public InequalityConstrainedSolver {
   public:
      BoxLPSolver() = default;
      ~BoxLPSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, Subproblem& subproblem, double trust_region_radius, const Vector<double>& initial_point,
         Direction& direction, const Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   protected:
      std::vector<double> variable_lower_bounds;
      std::vector<double> variable_upper_bounds;
      BoxLPSolverWorkspace evaluation_space{};
   };
} // namespace

#endif // UNO_BOXLPSOLVER_H