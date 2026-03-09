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

   // Lagrangian gradient ρ ∇f(x_k) - ∇c(x_k) y_k - z_k
   void Model::evaluate_lagrangian_gradient(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier,
         Evaluations& evaluations, Vector<double>& lagrangian_gradient) const {
      lagrangian_gradient.fill(0.);

      // - ∇c(x_k) λ_k
      // TODO test whether λ_k != 0
      evaluations.evaluate_jacobian(*this, primals);
      evaluations.compute_jacobian_transposed_vector_product(multipliers.constraints, lagrangian_gradient);
      lagrangian_gradient.scale(-1.);

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
      for (size_t variable_index: Range(this->number_variables)) {
         x[variable_index] = std::max(std::min(x[variable_index], this->variable_upper_bound(variable_index)), this->variable_lower_bound(variable_index));
      }
   }

   bool Model::is_constrained() const {
      return (0 < this->number_constraints);
   }

   // individual constraint violation
   double Model::constraint_violation(double constraint_value, size_t constraint_index) const {
      const double lower_bound_violation = std::max(0., this->constraint_lower_bound(constraint_index) - constraint_value);
      const double upper_bound_violation = std::max(0., constraint_value - this->constraint_upper_bound(constraint_index));
      return std::max(lower_bound_violation, upper_bound_violation);
   }

   void Model::find_fixed_variables(Vector<size_t>& fixed_variables) const {
      fixed_variables.reserve(this->number_variables);

      for (size_t variable_index: Range(this->number_variables)) {
         if (this->variable_lower_bound(variable_index) == this->variable_upper_bound(variable_index)) {
            WARNING << "Variable x" << variable_index << " has identical bounds\n";
            fixed_variables.emplace_back(variable_index);
         }
      }
   }

   void Model::partition_constraints(std::vector<size_t>& equality_constraints, std::vector<size_t>& inequality_constraints) const {
      equality_constraints.reserve(this->number_constraints);
      inequality_constraints.reserve(this->number_constraints);

      for (size_t constraint_index: Range(this->number_constraints)) {
         const double lower_bound = this->constraint_lower_bound(constraint_index);
         const double upper_bound = this->constraint_upper_bound(constraint_index);
         if (lower_bound == upper_bound) {
            equality_constraints.emplace_back(constraint_index);
         }
         else if (!is_finite(lower_bound) && !is_finite(upper_bound)) {
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