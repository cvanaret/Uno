// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "DoglegWorkspace.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/LinearSystem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   void DoglegWorkspace::initialize_memory(const Subproblem& subproblem) {
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

   void DoglegWorkspace::compute_newton_step(const Subproblem& subproblem, Direction& direction,
         DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      if (warmstart_information.new_iterate) {
         // g = ∇f(xₖ)
         subproblem.problem.evaluate_objective_gradient(subproblem.current_iterate, this->objective_gradient.data(),
            current_evaluations);
         // form matrix
         auto& linear_system = linear_solver.get_linear_system();
         Statistics statistics{};
         subproblem.evaluate_lagrangian_hessian(statistics, linear_system.matrix_values.data());
         std::cout << "MATRIX: " << linear_system.matrix_values << '\n';
         linear_solver.do_numerical_factorization(true);
         const Inertia inertia = linear_solver.get_inertia();
         std::cout << "INERTIA = " << inertia << '\n';
         // form RHS
         subproblem.problem.evaluate_objective_gradient(subproblem.current_iterate, linear_system.rhs.data(), current_evaluations);
         linear_system.rhs.scale(-1.);
         std::cout << "RHS: " << linear_system.rhs << '\n';
         linear_solver.solve_indefinite_system(linear_system.solution.data());
         std::cout << "DOGLEG: NEWTON SOLUTION: " << linear_system.solution << '\n';
         subproblem.assemble_primal_dual_direction(linear_system.solution, direction);
         this->newton_step = direction.primals;
         this->newton_step_squared_norm = dot(this->newton_step, this->newton_step);
      }
   }

   void DoglegWorkspace::compute_dogleg(const Subproblem& subproblem, Direction& /*direction*/, Evaluations& /*current_evaluations*/,
         const WarmstartInformation& warmstart_information) {
      this->compute_cauchy_step(subproblem, warmstart_information);
      // TODO
   }

   // private member functions

   void DoglegWorkspace::compute_cauchy_step(const Subproblem& subproblem, const WarmstartInformation& warmstart_information) {
      if (warmstart_information.new_iterate) {
         // gᵀ g
         this->objective_gradient_squared_norm = dot(this->objective_gradient, this->objective_gradient);
         // B g = Hₖ ∇f(xₖ)
         subproblem.compute_hessian_vector_product(subproblem.current_iterate.primals.data(), this->objective_gradient.data(),
            this->hessian_gradient_product.data());
         // gᵀ B g = ∇f(xₖ)ᵀ Hₖ ∇f(xₖ)
         this->hessian_directional_derivative = dot(this->objective_gradient, this->hessian_gradient_product);
         if (this->hessian_directional_derivative <= 0.) {
            throw std::runtime_error("The Hessian is not positive definite");
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