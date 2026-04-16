// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <utility>
#include "Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   // abstract Problem class
   Model::Model(std::string name, size_t number_variables, size_t number_constraints, double objective_sign,
                double lagrangian_sign_convention) :
         name(std::move(name)), number_variables(number_variables), number_constraints(number_constraints),
         optimization_sense(objective_sign), lagrangian_sign_convention(lagrangian_sign_convention) {
   }

   bool Model::has_inequality_constraints() const {
      const auto& constraints_lower_bounds = this->get_constraints_lower_bounds();
      const auto& constraints_upper_bounds = this->get_constraints_upper_bounds();
      for (size_t constraint_index: Range(this->number_constraints)) {
         if (constraints_lower_bounds[constraint_index] < constraints_upper_bounds[constraint_index]) {
            return true;
         }
      }
      return false;
   }

   bool Model::has_bound_constraints() const {
      const auto& variables_lower_bounds = this->get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->get_variables_upper_bounds();
      if (std::any_of(variables_lower_bounds.begin(), variables_lower_bounds.end(), is_finite<double>) ||
            std::any_of(variables_upper_bounds.begin(), variables_upper_bounds.end(), is_finite<double>)) {
         return true;
            }
      return false;
   }

   // Lagrangian gradient ρ ∇f(x_k) - ∇c(x_k) y_k - z_k
   void Model::evaluate_lagrangian_gradient(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier,
         Evaluations& evaluations, Vector<double>& lagrangian_gradient) const {
      lagrangian_gradient.fill(0.);

      // - ∇c(x_k) λ_k
      if (0 < this->number_constraints) {
         // TODO test whether λ_k != 0
         evaluations.evaluate_jacobian(*this, primals);
         evaluations.compute_jacobian_transposed_vector_product(*this, multipliers.constraints.data(), lagrangian_gradient.data());
         lagrangian_gradient.scale(-1.);
      }

      // ∇f(x_k)
      if (objective_multiplier != 0.) {
         evaluations.evaluate_objective_gradient(*this, primals);
         view(lagrangian_gradient, 0, this->number_variables) += objective_multiplier * evaluations.objective_gradient;
      }

      // z_k
      view(lagrangian_gradient, 0, this->number_variables) -= view(multipliers.lower_bounds, 0, this->number_variables);
      view(lagrangian_gradient, 0, this->number_variables) -= view(multipliers.upper_bounds, 0, this->number_variables);
   }

   void Model::project_onto_variable_bounds(Vector<double>& x) const {
      const auto& variables_lower_bounds = this->get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->get_variables_upper_bounds();
      for (size_t variable_index: Range(this->number_variables)) {
         x[variable_index] = std::max(std::min(x[variable_index], variables_upper_bounds[variable_index]),
            variables_lower_bounds[variable_index]);
      }
   }

   bool Model::is_constrained() const {
      return (0 < this->number_constraints);
   }

   // individual constraint violation
   double Model::constraint_violation(double constraint_value, size_t constraint_index) const {
      const double lower_bound_violation = std::max(0., this->get_constraints_lower_bounds()[constraint_index] - constraint_value);
      const double upper_bound_violation = std::max(0., constraint_value - this->get_constraints_upper_bounds()[constraint_index]);
      return std::max(lower_bound_violation, upper_bound_violation);
   }

   void Model::find_fixed_variables(Vector<size_t>& fixed_variables) const {
      fixed_variables.reserve(this->number_variables);

      const auto& variables_lower_bounds = this->get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->get_variables_upper_bounds();
      for (size_t variable_index: Range(this->number_variables)) {
         if (variables_lower_bounds[variable_index] == variables_upper_bounds[variable_index]) {
            WARNING << "Variable x" << variable_index << " has identical bounds\n";
            fixed_variables.emplace_back(variable_index);
         }
      }
   }

   void Model::partition_constraints(std::vector<size_t>& equality_constraints, std::vector<size_t>& inequality_constraints) const {
      equality_constraints.reserve(this->number_constraints);
      inequality_constraints.reserve(this->number_constraints);

      const auto& constraints_lower_bounds = this->get_constraints_lower_bounds();
      const auto& constraints_upper_bounds = this->get_constraints_upper_bounds();
      for (size_t constraint_index: Range(this->number_constraints)) {
         if (constraints_lower_bounds[constraint_index] == constraints_upper_bounds[constraint_index]) {
            equality_constraints.emplace_back(constraint_index);
         }
         else if (is_infinite(constraints_lower_bounds[constraint_index]) && is_infinite(constraints_upper_bounds[constraint_index])) {
            WARNING << "Constraint c" << constraint_index << " has no bounds\n";
            // count the constraint as inequality to avoid reindexing of the constraints
            inequality_constraints.emplace_back(constraint_index);
         }
         else {
            inequality_constraints.emplace_back(constraint_index);
         }
      }
   }
} // namespace