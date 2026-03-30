// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "DoglegWorkspace.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "symbolic/Sum.hpp"
#include "tools/Logger.hpp"

namespace uno {
   DoglegWorkspace::DoglegWorkspace(const Options& options):
         newton_solver(options) {
   }

   void DoglegWorkspace::initialize_memory(const Subproblem& subproblem) {
      // Newton solver
      this->newton_solver.initialize_memory(subproblem);
      this->initial_point.resize(subproblem.number_variables);
      // Newton step
      this->g.resize(subproblem.number_variables);
      this->newton_step.resize(subproblem.number_variables);
      // Cauchy step
      this->Hg.resize(subproblem.number_variables);
      this->cauchy_step.resize(subproblem.number_variables);
   }

   double DoglegWorkspace::compute_hessian_quadratic_form(const Subproblem& subproblem, const Vector<double>& vector) const {
      return this->newton_solver.get_workspace().compute_hessian_quadratic_form(subproblem, vector);
   }

   void DoglegWorkspace::compute_newton_step(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) {
      if (warmstart_information.new_iterate) {
         this->newton_solver.solve(statistics, subproblem, INF<double>, this->initial_point, direction, current_evaluations,
            warmstart_information);
         this->newton_step = direction.primals;
         this->newton_step_squared_norm = dot(this->newton_step, this->newton_step);
      }
      DEBUG << "Newton step: " << this->newton_step << '\n';
   }

   void DoglegWorkspace::compute_dogleg(const Subproblem& subproblem, double trust_region_radius, Direction& direction,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) {
      this->compute_cauchy_step(subproblem, current_evaluations, warmstart_information);
      DEBUG << "Cauchy step: " << this->cauchy_step << '\n';
      double squared_norm_cauchy = dot(this->cauchy_step, this->cauchy_step);
      if (trust_region_radius*trust_region_radius <= squared_norm_cauchy) {
         DEBUG << "The Cauchy step is outside the trust region. Returning a clipped direction\n";
         direction.primals = this->cauchy_step;
         direction.primals.scale(trust_region_radius/std::sqrt(squared_norm_cauchy));
         direction.norm = trust_region_radius;
         DEBUG << "Scaled Cauchy direction: " << direction.primals << '\n';
         return;
      }

      DEBUG << "Computing the dogleg step\n";
      // define temporary vectors
      Vector<double> u(subproblem.number_variables);
      u = this->newton_step - this->cauchy_step;
      Vector<double> v(subproblem.number_variables);
      v = 2.*this->cauchy_step - this->newton_step;
      const double a = dot(u, u);
      const double b = 2.*dot(u, v);
      const double c = dot(v, v) - trust_region_radius*trust_region_radius;
      DEBUG << "(a, b, c) = " << a << ", " << b << ", " << c << '\n';
      const double tau = DoglegWorkspace::compute_positive_root_quadratic_equation(a, b, c);
      DEBUG << "tau = " << tau << '\n';
      if (1. <= tau && tau <= 2.) {
         direction.primals = this->cauchy_step + (1. - tau)*(this->newton_step - this->cauchy_step);
         direction.norm = norm_2(view(direction.primals, 0, subproblem.problem.get_number_original_variables()));
      }
      else {
         throw std::runtime_error("The dogleg step could not be computed");
      }
   }

   // private member functions

   void DoglegWorkspace::compute_cauchy_step(const Subproblem& subproblem, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      if (warmstart_information.new_iterate) {
         subproblem.problem.evaluate_objective_gradient(subproblem.current_iterate, this->g.data(), current_evaluations);
         DEBUG << "g = " << this->g << '\n';
         // gᵀ g
         this->g_squared_norm = dot(this->g, this->g);
         // H g
         subproblem.compute_hessian_vector_product(subproblem.current_iterate.primals.data(), this->g.data(),
            this->Hg.data());
         DEBUG << "H g = " << this->Hg << '\n';
         // gᵀ H g
         this->gHg = dot(this->g, this->Hg);
         DEBUG << "g^T H g = " << this->gHg << '\n';
         if (this->gHg <= 0.) {
            throw std::runtime_error("The Hessian is not positive definite");
         }
         // Cauchy step: d_C = - (gᵀ g)/(gᵀ B g) g
         this->cauchy_step = this->g;
         const double scaling_factor = -this->g_squared_norm / this->gHg;
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