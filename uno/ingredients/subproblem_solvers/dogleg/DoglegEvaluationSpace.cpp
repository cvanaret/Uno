// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "DoglegEvaluationSpace.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolver.hpp"
#include "optimization/Direction.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"

namespace uno {
   void DoglegEvaluationSpace::initialize_memory(const Subproblem& subproblem) {
      // Newton step
      this->objective_gradient.resize(subproblem.number_variables);
      this->newton_step.resize(subproblem.number_variables);
      // Cauchy step
      this->hessian_gradient_product.resize(subproblem.number_variables);
      this->cauchy_step.resize(subproblem.number_variables);
   }

   void DoglegEvaluationSpace::evaluate_constraint_jacobian(const OptimizationProblem& /*problem*/, Iterate& /*iterate*/) {
      // do nothing
   }

   void DoglegEvaluationSpace::compute_constraint_jacobian_vector_product(const Vector<double>& /*vector*/, Vector<double>& result) const {
      result.fill(0.);
   }

   void DoglegEvaluationSpace::compute_constraint_jacobian_transposed_vector_product(const Vector<double>& /*vector*/,
         Vector<double>& result) const {
      result.fill(0.);
   }

   double DoglegEvaluationSpace::compute_hessian_quadratic_product(const Subproblem& /*subproblem*/,
         const Vector<double>& /*vector*/) const {
      throw std::runtime_error("Not implemented yet");
   }

   void DoglegEvaluationSpace::evaluate_objective_gradient(const Subproblem& subproblem, const WarmstartInformation& warmstart_information) {
      if (warmstart_information.objective_changed) {
         subproblem.problem.evaluate_objective_gradient(subproblem.current_iterate, this->objective_gradient.data());
      }
   }

   void DoglegEvaluationSpace::compute_newton_step(Statistics& statistics, const Subproblem& subproblem,
         SymmetricIndefiniteLinearSolver<double>& linear_solver, Direction& direction, const WarmstartInformation& warmstart_information) {
      if (warmstart_information.objective_changed) {
         // g = ∇f(x_k)
         subproblem.problem.evaluate_objective_gradient(subproblem.current_iterate, this->objective_gradient.data());
         linear_solver.solve_indefinite_system(statistics, subproblem, direction, warmstart_information);
         this->newton_step = direction.primals;
         this->newton_step_squared_norm = dot(this->newton_step, this->newton_step);
      }
   }

   void DoglegEvaluationSpace::compute_dogleg(const Subproblem& subproblem, Direction& /*direction*/,
         const WarmstartInformation& warmstart_information) {
      this->compute_cauchy_step(subproblem, warmstart_information);
      // TODO
   }

   // private member functions

   void DoglegEvaluationSpace::compute_cauchy_step(const Subproblem& subproblem, const WarmstartInformation& warmstart_information) {
      if (warmstart_information.objective_changed) {
         // g^T g
         this->objective_gradient_squared_norm = dot(this->objective_gradient, this->objective_gradient);
         // B g = H_k ∇f(x_k)
         subproblem.compute_hessian_vector_product(subproblem.current_iterate.primals.data(), this->objective_gradient.data(),
            this->hessian_gradient_product.data());
         // g^T B g = ∇f(x_k)^T H_k ∇f(x_k)
         this->hessian_quadratic_product = dot(this->objective_gradient, this->hessian_gradient_product);
         if (this->hessian_quadratic_product <= 0.) {
            throw std::runtime_error("The objective Hessian is not positive definite");
         }
         // Cauchy step: d_C = - (g^T g)/(g^T B g) g
         this->cauchy_step = this->objective_gradient;
         const double scaling_factor = -this->objective_gradient_squared_norm / this->hessian_quadratic_product;
         this->cauchy_step.scale(scaling_factor);
      }
   }

   // find the positive real root to ax^2 + bx + c = 0
   double DoglegEvaluationSpace::compute_positive_root_quadratic_equation(double a, double b, double c) {
      // avoid catastrophic cancellation
      const double delta = b*b - 4.*a*c;
      if (delta < 0.) {
         throw std::runtime_error("No real root");
      }
      // TODO check denominator
      return (2.*c)/(-b - std::sqrt(delta));
   }
} // namespace