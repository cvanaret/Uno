// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <utility>
#include "Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   // abstract Problem class
   Model::Model(std::string name, size_t number_variables, size_t number_constraints, double objective_sign) :
         name(std::move(name)), number_variables(number_variables), number_constraints(number_constraints),
         optimization_sense(objective_sign) {
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