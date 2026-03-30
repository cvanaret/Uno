// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DOGLEGWORKSPACE_H
#define UNO_DOGLEGWORKSPACE_H

#include "../SolverWorkspace.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   // forward declaration
   class Direction;
   class Evaluations;
   class Statistics;
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class WarmstartInformation;

   class DoglegWorkspace: public SolverWorkspace {
   public:
      DoglegWorkspace() = default;
      ~DoglegWorkspace() override = default;

      // Newton step
      Vector<double> objective_gradient{};
      Vector<double> newton_step{};
      double newton_step_squared_norm{INF<double>};
      // Cauchy step
      double objective_gradient_squared_norm{INF<double>};
      Vector<double> hessian_gradient_product{};
      double hessian_directional_derivative{INF<double>};
      Vector<double> cauchy_step{};

      void initialize_memory(const Subproblem& subproblem);

      [[nodiscard]] double compute_hessian_quadratic_form(const Subproblem& subproblem, const Vector<double>& vector) const override;

      void compute_newton_step(const Subproblem& subproblem, Direction& direction, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information);
      void compute_dogleg(const Subproblem& subproblem, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information);

   private:
      void compute_cauchy_step(const Subproblem& subproblem, const WarmstartInformation& warmstart_information);
      [[nodiscard]] static double compute_positive_root_quadratic_equation(double a, double b, double c);
   };
} // namespace

#endif // UNO_DOGLEGWORKSPACE_H