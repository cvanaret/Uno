// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HiGHSQuadraticProgram.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   void HiGHSQuadraticProgram::initialize_memory(const Subproblem& subproblem) {
      this->number_variables = subproblem.number_variables;
      this->number_constraints = subproblem.number_constraints;
      this->number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      this->workspace.initialize_memory(subproblem);
   }

   void HiGHSQuadraticProgram::fill(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
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

   void HiGHSQuadraticProgram::fill(const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values,
         const Vector<uno_int>& hessian_row_indices, const Vector<uno_int>& hessian_column_indices,
         const Vector<double>& hessian_values,
         const std::vector<double>& variables_lower_bounds, const std::vector<double>& variables_upper_bounds,
         const std::vector<double>& constraints_lower_bounds, const std::vector<double>& constraints_upper_bounds) {
      // infer the dimensions from the data
      this->number_variables = linear_objective.size();
      this->number_constraints = constraints_lower_bounds.size();

      // allocate the HighsModel and convert the COO Jacobian/Hessian to HiGHS' CSC layout
      this->workspace.set_from_coo(this->number_variables, this->number_constraints, linear_objective,
         jacobian_row_indices, jacobian_column_indices, jacobian_values,
         hessian_row_indices, hessian_column_indices, hessian_values);

      // variable and constraint bounds (HiGHS uses its own infinity, so no clamping is required)
      for (size_t variable_index: Range(this->number_variables)) {
         this->workspace.model.lp_.col_lower_[variable_index] = variables_lower_bounds[variable_index];
         this->workspace.model.lp_.col_upper_[variable_index] = variables_upper_bounds[variable_index];
      }
      for (size_t constraint_index: Range(this->number_constraints)) {
         this->workspace.model.lp_.row_lower_[constraint_index] = constraints_lower_bounds[constraint_index];
         this->workspace.model.lp_.row_upper_[constraint_index] = constraints_upper_bounds[constraint_index];
      }
   }

   SolverWorkspace& HiGHSQuadraticProgram::get_workspace() {
      return this->workspace;
   }
} // namespace
