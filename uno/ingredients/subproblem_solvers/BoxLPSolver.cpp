// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "BoxLPSolver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"

namespace uno {
   void BoxLPSolver::initialize_memory(const Subproblem& subproblem) {
      this->variable_lower_bounds.resize(subproblem.number_variables);
      this->variable_upper_bounds.resize(subproblem.number_variables);
      this->evaluation_space.objective_gradient.resize(subproblem.number_variables);
   }

   void BoxLPSolver::solve(Statistics& /*statistics*/, Subproblem& subproblem, const Vector<double>& /*initial_point*/,
         Direction& direction, const WarmstartInformation& /*warmstart_information*/) {
      if (0 < subproblem.number_constraints) {
         throw std::runtime_error("BoxLPSolver cannot solve problems with general constraints");
      }
      // compute the objective gradient
      subproblem.problem.evaluate_objective_gradient(subproblem.current_iterate, this->evaluation_space.objective_gradient.data());

      // compute the variables bounds
      subproblem.set_variables_bounds(this->variable_lower_bounds, this->variable_upper_bounds);

      // move the variables to one of their bounds
      direction.subproblem_objective = 0.;
      for (size_t variable_index: Range(subproblem.number_variables)) {
         if (0. < this->evaluation_space.objective_gradient[variable_index]) {
            direction.primals[variable_index] = this->variable_lower_bounds[variable_index];
            direction.multipliers.lower_bounds[variable_index] = this->evaluation_space.objective_gradient[variable_index];
            if (!is_finite(this->variable_lower_bounds[variable_index])) {
               direction.status = SubproblemStatus::UNBOUNDED_PROBLEM;
            }
         }
         else if (this->evaluation_space.objective_gradient[variable_index] < 0.) {
            direction.primals[variable_index] = this->variable_upper_bounds[variable_index];
            direction.multipliers.upper_bounds[variable_index] = this->evaluation_space.objective_gradient[variable_index];
            if (!is_finite(this->variable_upper_bounds[variable_index])) {
               direction.status = SubproblemStatus::UNBOUNDED_PROBLEM;
            }
         }
         else {
            direction.primals[variable_index] = 0.;
            direction.multipliers.lower_bounds[variable_index] = direction.multipliers.upper_bounds[variable_index] = 0.;
         }
         direction.subproblem_objective += this->evaluation_space.objective_gradient[variable_index] * direction.primals[variable_index];
      }
   }

   EvaluationSpace& BoxLPSolver::get_evaluation_space() {
      return this->evaluation_space;
   }
} // namespace