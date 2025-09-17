// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BOXLPSOLVER_H
#define UNO_BOXLPSOLVER_H

#include <vector>
#include "InequalityConstrainedSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/EvaluationSpace.hpp"

namespace uno {
   class BoxLPSolverEvaluationSpace: public EvaluationSpace {
   public:
      BoxLPSolverEvaluationSpace() = default;

      void evaluate_constraint_jacobian(const OptimizationProblem& /*problem*/, Iterate& /*iterate*/) override { }
      void compute_constraint_jacobian_vector_product(const Vector<double>& /*vector*/, Vector<double>& /*result*/) const override { }
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& /*vector*/,
         Vector<double>& /*result*/) const override { }
      [[nodiscard]] double compute_hessian_quadratic_product(const Vector<double>& /*vector*/) const override {
         return 0.;
      }

      Vector<double> objective_gradient;
   };

   class BoxLPSolver: public InequalityConstrainedSolver {
   public:
      BoxLPSolver() = default;
      ~BoxLPSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, Subproblem& subproblem, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] EvaluationSpace& get_evaluation_space() override;

   protected:
      std::vector<double> variable_lower_bounds;
      std::vector<double> variable_upper_bounds;
      BoxLPSolverEvaluationSpace evaluation_space{};
   };
} // namespace

#endif // UNO_BOXLPSOLVER_H