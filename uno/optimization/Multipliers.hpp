// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MULTIPLIERS_H
#define UNO_MULTIPLIERS_H

#include "linear_algebra/Vector.hpp"

struct Multipliers {
   std::vector<double> lower_bounds{}; /*!< Multipliers of the lower bound constraints */
   std::vector<double> upper_bounds{}; /*!< Multipliers of the lower bound constraints */
   std::vector<double> constraints{}; /*!< Multipliers of the general constraints */
   double objective{1.};

   Multipliers(size_t number_variables, size_t number_constraints);
   [[nodiscard]] bool not_all_zero(size_t number_variables, double tolerance) const;
};

inline Multipliers::Multipliers(size_t number_variables, size_t number_constraints) : lower_bounds(number_variables),
      upper_bounds(number_variables), constraints(number_constraints) {
}

inline bool Multipliers::not_all_zero(size_t number_variables, double tolerance) const {
   // constraint multipliers
   for (double multiplier_j: this->constraints) {
      if (tolerance < std::abs(multiplier_j)) {
         return true;
      }
   }
   // bound multipliers
   for (size_t i: Range(number_variables)) {
      if (tolerance < std::abs(this->lower_bounds[i] + this->upper_bounds[i])) {
         return true;
      }
   }
   return false;
}

#endif // UNO_MULTIPLIERS_H
