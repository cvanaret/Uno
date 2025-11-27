// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "DoglegMethod.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "options/Options.hpp"

namespace uno {
   DoglegMethod::DoglegMethod(const Options& options):
         linear_solver_name(options.get_string("linear_solver")) {
   }

   void DoglegMethod::initialize_memory(const Subproblem& subproblem) {
      this->linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->linear_solver_name);
      this->linear_solver->initialize_hessian(subproblem);

      this->evaluation_space.initialize_memory(subproblem);
   }

   void DoglegMethod::solve(Statistics& statistics, Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& /*initial_point*/, Direction& direction, const WarmstartInformation& warmstart_information) {
      this->evaluation_space.evaluate_objective_gradient(subproblem, warmstart_information);

      const double squared_trust_region_radius = std::pow(trust_region_radius, 2.);
      // first try the Newton step. This is the solution if within the trust region
      this->evaluation_space.compute_newton_step(statistics, subproblem, *this->linear_solver, direction, warmstart_information);
      if (this->evaluation_space.newton_step_squared_norm <= squared_trust_region_radius) {
         return;
      }
      // if the trust region constraint is violated, compute the dogleg path: the broken path between the Cauchy step
      // and the Newton step
      this->evaluation_space.compute_dogleg(subproblem, direction, warmstart_information);
   }

   [[nodiscard]] EvaluationSpace& DoglegMethod::get_evaluation_space() {
      return this->evaluation_space;
   }
} // namespace