// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <iostream>
#include <utility>
#include "Model.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // abstract Problem class
   Model::Model(std::string name, size_t number_variables, size_t number_constraints, double objective_sign) :
         name(std::move(name)), number_variables(number_variables), number_constraints(number_constraints), objective_sign(objective_sign) {
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
} // namespace
