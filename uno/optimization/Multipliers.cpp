// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Multipliers.hpp"

namespace uno {
   Multipliers::Multipliers(size_t number_variables, size_t number_constraints) : lower_bounds(number_variables),
         upper_bounds(number_variables), constraints(number_constraints) {
   }

   void Multipliers::reset() {
      this->constraints.fill(0.);
      this->lower_bounds.fill(0.);
      this->upper_bounds.fill(0.);
   }

   bool Multipliers::not_all_zero(size_t number_variables, double tolerance) const {
      // constraint multipliers
      for (double multiplier_j: this->constraints) {
         if (tolerance < std::abs(multiplier_j)) {
            return true;
         }
      }
      // bound multipliers
      for (size_t variable_index: Range(number_variables)) {
         if (tolerance < std::abs(this->lower_bounds[variable_index] + this->upper_bounds[variable_index])) {
            return true;
         }
      }
      return false;
   }
} // namespace