#include <vector>
#include "Subproblem.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"

namespace uno {
   Subproblem::Subproblem(const OptimizationProblem &problem, Iterate &current_iterate, const Multipliers &current_multipliers,
      HessianModel &hessian_model, RegularizationStrategy<double> &regularization_strategy, double trust_region_radius):
         number_variables(problem.number_variables), number_constraints(problem.number_constraints),
         problem(problem), current_iterate(current_iterate), current_multipliers(current_multipliers), hessian_model(hessian_model),
         regularization_strategy(regularization_strategy), trust_region_radius(trust_region_radius) {
   }

   void Subproblem::evaluate_functions(SparseVector<double>& linear_objective, std::vector<double>& constraints,
         RectangularMatrix<double>& constraint_jacobian, const WarmstartInformation& warmstart_information) const {
      if (warmstart_information.objective_changed) {
         this->problem.evaluate_objective_gradient(this->current_iterate, linear_objective);
      }
      if (warmstart_information.constraints_changed) {
         this->problem.evaluate_constraints(this->current_iterate, constraints);
         this->problem.evaluate_constraint_jacobian(this->current_iterate, constraint_jacobian);
      }
   }

   void Subproblem::compute_regularized_hessian(Statistics& statistics, SymmetricMatrix<size_t, double>& hessian,
         const WarmstartInformation& warmstart_information) const {
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         // evaluate the Lagrangian Hessian of the problem at the current primal-dual point
         this->problem.evaluate_lagrangian_hessian(statistics, this->hessian_model, this->current_iterate.primals,
            this->current_multipliers, hessian);
         // regularize the Hessian
         const Inertia expected_inertia{this->problem.get_number_original_variables(), 0,
            this->problem.number_variables - this->problem.get_number_original_variables()};
         this->regularization_strategy.regularize_hessian(statistics, hessian, expected_inertia);
      }
   }

   void Subproblem::assemble_augmented_matrix(Statistics& statistics, SymmetricMatrix<size_t, double>& augmented_matrix,
         RectangularMatrix<double>& constraint_jacobian) const {
      // evaluate the Lagrangian Hessian of the problem at the current primal-dual point
      this->problem.evaluate_lagrangian_hessian(statistics, this->hessian_model, this->current_iterate.primals,
         this->current_multipliers, augmented_matrix);
      augmented_matrix.set_dimension(this->problem.number_variables + this->problem.number_constraints);

      // Jacobian of general constraints
      for (size_t column_index: Range(this->problem.number_constraints)) {
         for (const auto [row_index, derivative]: constraint_jacobian[column_index]) {
            augmented_matrix.insert(derivative, row_index, this->problem.number_variables + column_index);
         }
         augmented_matrix.finalize_column(column_index);
      }
   }

   void Subproblem::assemble_augmented_rhs(const SparseVector<double>& objective_gradient, const std::vector<double>& constraints,
         RectangularMatrix<double>& constraint_jacobian, Vector<double>& rhs) const {
      rhs.fill(0.);

      // objective gradient
      for (const auto [variable_index, derivative]: objective_gradient) {
         rhs[variable_index] -= derivative;
      }

      // constraint: evaluations and gradients
      for (size_t constraint_index: Range(number_constraints)) {
         // Lagrangian
         if (current_multipliers.constraints[constraint_index] != 0.) {
            for (const auto [variable_index, derivative]: constraint_jacobian[constraint_index]) {
               rhs[variable_index] += current_multipliers.constraints[constraint_index] * derivative;
            }
         }
         // constraints
         rhs[number_variables + constraint_index] = -constraints[constraint_index];
      }
      DEBUG2 << "RHS: "; print_vector(DEBUG2, view(rhs, 0, number_variables + number_constraints)); DEBUG << '\n';
   }

   void Subproblem::set_variables_bounds(std::vector<double>& variables_lower_bounds, std::vector<double>& variables_upper_bounds,
         const WarmstartInformation& warmstart_information) const {
      if (warmstart_information.variable_bounds_changed) {
         // bounds of original variables intersected with trust region
         for (size_t variable_index: Range(this->problem.get_number_original_variables())) {
            variables_lower_bounds[variable_index] = std::max(-this->trust_region_radius,
               this->problem.variable_lower_bound(variable_index) - this->current_iterate.primals[variable_index]);
            variables_upper_bounds[variable_index] = std::min(this->trust_region_radius,
               this->problem.variable_upper_bound(variable_index) - this->current_iterate.primals[variable_index]);
         }
         // bounds of additional variables (no trust region!)
         for (size_t variable_index: Range(this->problem.get_number_original_variables(), this->problem.number_variables)) {
            variables_lower_bounds[variable_index] = this->problem.variable_lower_bound(variable_index) - this->current_iterate.primals[variable_index];
            variables_upper_bounds[variable_index] = this->problem.variable_upper_bound(variable_index) - this->current_iterate.primals[variable_index];
         }
      }
   }
} // namespace