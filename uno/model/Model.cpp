// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <utility>
#include "Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/LagrangianGradient.hpp"
#include "optimization/Multipliers.hpp"

namespace uno {
   Model::Model(std::string name, size_t number_variables, size_t number_constraints, double objective_sign) :
         name(std::move(name)), number_variables(number_variables), number_constraints(number_constraints), objective_sign(objective_sign) {
   }

   void Model::project_onto_variable_bounds(Vector<double>& x) const {
      for (size_t variable_index: Range(this->number_variables)) {
         x[variable_index] = std::max(std::min(x[variable_index], this->variable_upper_bound(variable_index)),
            this->variable_lower_bound(variable_index));
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

   // Lagrangian gradient split in two parts: objective contribution and constraints' contribution
   void Model::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, const Iterate& iterate,
         const Multipliers& multipliers) const {
      lagrangian_gradient.objective_contribution.fill(0.);
      lagrangian_gradient.constraints_contribution.fill(0.);

      // objective gradient
      /*
      for (auto [variable_index, derivative]: iterate.evaluations.objective_gradient) {
         lagrangian_gradient.objective_contribution[variable_index] += derivative;
      }
      */

      // constraints
      for (size_t constraint_index: Range(this->number_constraints)) {
         if (multipliers.constraints[constraint_index] != 0.) {
            for (auto [variable_index, derivative]: iterate.evaluations.constraint_jacobian[constraint_index]) {
               lagrangian_gradient.constraints_contribution[variable_index] -= multipliers.constraints[constraint_index] * derivative;
            }
         }
      }

      // bound constraints of original variables
      for (size_t variable_index: Range(this->number_variables)) {
         lagrangian_gradient.constraints_contribution[variable_index] -= (multipliers.lower_bounds[variable_index] +
            multipliers.upper_bounds[variable_index]);
      }
   }
} // namespace
