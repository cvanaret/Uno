// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MULTIPLIERS_H
#define UNO_MULTIPLIERS_H

#include "linear_algebra/Vector.hpp"

struct Multipliers {
   Vector<double> lower_bounds{}; /*!< Multipliers of the lower bound constraints */
   Vector<double> upper_bounds{}; /*!< Multipliers of the lower bound constraints */
   Vector<double> constraints{}; /*!< Multipliers of the general constraints */

   Multipliers(size_t number_variables, size_t number_constraints);
   Multipliers(Multipliers&& other) noexcept = default;
   Multipliers& operator=(const Multipliers& other) noexcept = default;
   Multipliers& operator=(Multipliers&& other) noexcept = default;

   void reset();
   [[nodiscard]] bool not_all_zero(size_t number_variables, double tolerance) const;
};

inline Multipliers::Multipliers(size_t number_variables, size_t number_constraints) : lower_bounds(number_variables),
      upper_bounds(number_variables), constraints(number_constraints) {
}

inline void Multipliers::reset() {
   this->constraints.fill(0.);
   this->lower_bounds.fill(0.);
   this->upper_bounds.fill(0.);
}

inline bool Multipliers::not_all_zero(size_t number_variables, double tolerance) const {
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

#endif // UNO_MULTIPLIERS_H
