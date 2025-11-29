// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DOGLEGEVALUATIONSPACE_H
#define UNO_DOGLEGEVALUATIONSPACE_H

#include "optimization/EvaluationSpace.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   // forward declaration
   class Direction;
   template <typename ElementType>
   class SymmetricIndefiniteLinearSolver;

   class DoglegEvaluationSpace: public EvaluationSpace {
   public:
      DoglegEvaluationSpace() = default;
      ~DoglegEvaluationSpace() override = default;

      // Newton step
      Vector<double> objective_gradient{};
      Vector<double> newton_step{};
      double newton_step_squared_norm{INF<double>};
      // Cauchy step
      double objective_gradient_squared_norm{INF<double>};
      Vector<double> hessian_gradient_product{};
      double hessian_quadratic_product{INF<double>};
      Vector<double> cauchy_step{};

      void initialize_memory(const Subproblem& subproblem);

      void evaluate_constraint_jacobian(const OptimizationProblem& /*problem*/, Iterate& /*iterate*/) override;
      void compute_constraint_jacobian_vector_product(const Vector<double>& /*vector*/, Vector<double>& result) const override;
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      [[nodiscard]] double compute_hessian_quadratic_product(const Subproblem& subproblem,
         const Vector<double>& vector) const override;

      void evaluate_objective_gradient(const Subproblem& subproblem, const WarmstartInformation& warmstart_information);
      void compute_newton_step(Statistics& statistics, const Subproblem& subproblem,
         SymmetricIndefiniteLinearSolver<double>& linear_solver, Direction& direction, const WarmstartInformation& warmstart_information);
      void compute_dogleg(const Subproblem& subproblem, Direction& direction, const WarmstartInformation& warmstart_information);

   private:
      void compute_cauchy_step(const Subproblem& subproblem, const WarmstartInformation& warmstart_information);
      [[nodiscard]] static double compute_positive_root_quadratic_equation(double a, double b, double c);
   };
} // namespace

#endif // UNO_DOGLEGEVALUATIONSPACE_H