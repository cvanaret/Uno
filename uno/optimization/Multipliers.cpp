// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "Multipliers.hpp"

namespace uno {
   Multipliers::Multipliers(size_t number_variables, size_t number_constraints) : lower_bounds(number_variables),
         upper_bounds(number_variables), constraints(number_constraints) {
   }

   void Multipliers::resize(size_t number_variables, size_t number_constraints) {
      this->constraints.resize(number_constraints);
      this->lower_bounds.resize(number_variables);
      this->upper_bounds.resize(number_variables);
   }

   void Multipliers::reset() {
      this->constraints.fill(0.);
      this->lower_bounds.fill(0.);
      this->upper_bounds.fill(0.);
   }
} // namespace