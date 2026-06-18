// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HiGHSQuadraticProgram.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/WarmstartInformation.hpp"

namespace uno {
   HiGHSQuadraticProgram::HiGHSQuadraticProgram(size_t number_variables, size_t number_constraints):
      QuadraticProgram(number_variables, number_constraints) {
   }

   void HiGHSQuadraticProgram::initialize_memory(const Subproblem& subproblem) {
      this->workspace.initialize_memory(subproblem);
   }

   void HiGHSQuadraticProgram::build(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) {
      // evaluate the functions and derivatives into the HiGHS model
      this->workspace.evaluate_functions(statistics, subproblem, current_evaluations, warmstart_information);

      // variable bounds
      if (warmstart_information.trust_region_changed) {
         subproblem.set_variables_bounds(this->workspace.model.lp_.col_lower_, this->workspace.model.lp_.col_upper_,
            trust_region_radius);
      }
      // constraint bounds
      if (warmstart_information.constraint_bounds_changed || warmstart_information.new_iterate) {
         subproblem.set_constraints_bounds(this->workspace.model.lp_.row_lower_, this->workspace.model.lp_.row_upper_,
            this->workspace.constraints);
      }
   }

   SolverWorkspace& HiGHSQuadraticProgram::get_workspace() {
      return this->workspace;
   }
} // namespace
