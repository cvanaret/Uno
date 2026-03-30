// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "DoglegMethod.hpp"
#include "ingredients/subproblem/Subproblem.hpp"

namespace uno {
   DoglegMethod::DoglegMethod(const Options& options):
         workspace(options) {
   }

   void DoglegMethod::initialize_memory(const Subproblem& subproblem) {
      this->workspace.initialize_memory(subproblem);
   }

   void DoglegMethod::solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& /*initial_point*/, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      const double squared_trust_region_radius = std::pow(trust_region_radius, 2.);
      // first try the Newton step. This is the solution if within the trust region
      this->workspace.compute_newton_step(statistics, subproblem, direction, current_evaluations, warmstart_information);
      if (this->workspace.newton_step_squared_norm <= squared_trust_region_radius) {
         return;
      }
      // if the trust region constraint is violated, compute the dogleg path: the broken path between the Cauchy step
      // and the Newton step
      this->workspace.compute_dogleg(subproblem, direction, current_evaluations, warmstart_information);
   }

   [[nodiscard]] SolverWorkspace& DoglegMethod::get_workspace() {
      return this->workspace;
   }
} // namespace