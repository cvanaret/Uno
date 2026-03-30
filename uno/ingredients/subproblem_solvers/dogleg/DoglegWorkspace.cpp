// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "DoglegWorkspace.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"

namespace uno {
   DoglegWorkspace::DoglegWorkspace(const Options& options):
         newton_solver(options) {
   }

   void DoglegWorkspace::initialize_memory(const Subproblem& subproblem) {
      // Newton solver
      this->newton_solver.initialize_memory(subproblem);
      // Newton step
      this->objective_gradient.resize(subproblem.number_variables);
      this->newton_step.resize(subproblem.number_variables);
      // Cauchy step
      this->hessian_gradient_product.resize(subproblem.number_variables);
      this->cauchy_step.resize(subproblem.number_variables);
   }

   double DoglegWorkspace::compute_hessian_quadratic_form(const Subproblem& /*subproblem*/, const Vector<double>& /*vector*/) const {
      throw std::runtime_error("Not implemented yet");
   }

   void DoglegWorkspace::compute_newton_step(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) {
      Vector<double> initial_point(0);
      this->newton_solver.solve(statistics, subproblem, INF<double>, initial_point, direction, current_evaluations,
         warmstart_information);
      this->newton_step = direction.primals;
      this->newton_step_squared_norm = dot(this->newton_step, this->newton_step);
   }

   void DoglegWorkspace::compute_dogleg(const Subproblem& subproblem, Direction& /*direction*/, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      this->compute_cauchy_step(subproblem, current_evaluations, warmstart_information);
      // TODO
   }

   // private member functions

   void DoglegWorkspace::compute_cauchy_step(const Subproblem& subproblem, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      if (warmstart_information.new_iterate) {
         subproblem.problem.evaluate_objective_gradient(subproblem.current_iterate, this->objective_gradient.data(),
            current_evaluations);
         std::cout << "∇f(xₖ) = " << this->objective_gradient << '\n';
         // gᵀ g
         this->objective_gradient_squared_norm = dot(this->objective_gradient, this->objective_gradient);
         // B g = Hₖ ∇f(xₖ)
         subproblem.compute_hessian_vector_product(subproblem.current_iterate.primals.data(), this->objective_gradient.data(),
            this->hessian_gradient_product.data());
         std::cout << "Hₖ ∇f(xₖ) = " << this->hessian_gradient_product << '\n';
         // gᵀ B g = ∇f(xₖ)ᵀ Hₖ ∇f(xₖ)
         this->hessian_directional_derivative = dot(this->objective_gradient, this->hessian_gradient_product);
         std::cout << "∇f(xₖ)ᵀ Hₖ ∇f(xₖ) = " << this->hessian_directional_derivative << '\n';
         if (this->hessian_directional_derivative < 0.) {
            throw std::runtime_error("The Hessian is not positive semi-definite");
         }
         // Cauchy step: d_C = - (gᵀ g)/(gᵀ B g) g
         this->cauchy_step = this->objective_gradient;
         const double scaling_factor = -this->objective_gradient_squared_norm / this->hessian_directional_derivative;
         this->cauchy_step.scale(scaling_factor);
      }
   }

   // find the positive real root to ax^2 + bx + c = 0
   double DoglegWorkspace::compute_positive_root_quadratic_equation(double a, double b, double c) {
      // avoid catastrophic cancellation
      const double delta = b*b - 4.*a*c;
      if (delta < 0.) {
         throw std::runtime_error("No real root");
      }
      // TODO check denominator
      return (2.*c)/(-b - std::sqrt(delta));
   }
} // namespace