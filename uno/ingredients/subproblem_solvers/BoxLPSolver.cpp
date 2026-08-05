// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "BoxLPSolver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "tools/Logger.hpp"

namespace uno {
   void BoxLPSolver::initialize_memory(const Subproblem& subproblem) {
      this->direction = Direction(subproblem.number_variables, subproblem.number_constraints);
      this->variable_lower_bounds.resize(subproblem.number_variables);
      this->variable_upper_bounds.resize(subproblem.number_variables);
      this->workspace.objective_gradient.resize(subproblem.number_variables);
   }

   Direction& BoxLPSolver::solve(Statistics& /*statistics*/, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& /*initial_point*/, Evaluations& current_evaluations,
         const WarmstartInformation& /*warmstart_information*/) {
      if (0 < subproblem.number_constraints) {
         throw std::runtime_error("BoxLPSolver cannot solve problems with general constraints");
      }
      this->direction.reset();
      // compute the objective gradient
      this->workspace.objective_gradient.fill(0.);
      subproblem.problem.evaluate_objective_gradient(subproblem.current_iterate, this->workspace.objective_gradient.data(),
        current_evaluations);

      // compute the variables bounds
      subproblem.set_variables_bounds(this->variable_lower_bounds, this->variable_upper_bounds, trust_region_radius);

      // move the variables to one of their bounds
      this->direction.subproblem_objective = 0.;
      for (size_t variable_index: Range(subproblem.number_variables)) {
         if (0. < this->workspace.objective_gradient[variable_index]) {
            this->direction.primals[variable_index] = this->variable_lower_bounds[variable_index];
            this->direction.multipliers.lower_bounds[variable_index] = this->workspace.objective_gradient[variable_index];
            if (is_infinite(this->variable_lower_bounds[variable_index])) {
               this->direction.status = SubproblemStatus::UNBOUNDED_PROBLEM;
            }
         }
         else if (this->workspace.objective_gradient[variable_index] < 0.) {
            this->direction.primals[variable_index] = this->variable_upper_bounds[variable_index];
            this->direction.multipliers.upper_bounds[variable_index] = this->workspace.objective_gradient[variable_index];
            if (is_infinite(this->variable_upper_bounds[variable_index])) {
               this->direction.status = SubproblemStatus::UNBOUNDED_PROBLEM;
            }
         }
         else {
            this->direction.primals[variable_index] = 0.;
            this->direction.multipliers.lower_bounds[variable_index] = this->direction.multipliers.upper_bounds[variable_index] = 0.;
         }
         this->direction.subproblem_objective += this->workspace.objective_gradient[variable_index] * this->direction.primals[variable_index];
      }
      return this->direction;
   }

   bool BoxLPSolver::has_second_order_corrections() const {
      return false;
   }

   const Direction& BoxLPSolver::compute_second_order_correction(const Subproblem& /*subproblem*/, const Vector<double>& /*constraints_SOC*/) {
      throw std::runtime_error("No SOC implemented in BoxLPSolver");
   }

   SolverWorkspace& BoxLPSolver::get_workspace() {
      return this->workspace;
   }
} // namespace