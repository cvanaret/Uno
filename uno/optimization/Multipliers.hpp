// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MULTIPLIERS_H
#define UNO_MULTIPLIERS_H

#include <cmath>
#include "linear_algebra/Vector.hpp"

namespace uno {
   struct Multipliers {
      Vector<double> lower_bounds{}; /*!< Multipliers of the lower bound constraints */
      Vector<double> upper_bounds{}; /*!< Multipliers of the lower bound constraints */
      Vector<double> constraints{}; /*!< Multipliers of the general constraints */

      Multipliers(size_t number_variables, size_t number_constraints);
      Multipliers(const Multipliers& other) = default;
      Multipliers(Multipliers&& other) = default;
      Multipliers& operator=(const Multipliers& other) = default;
      Multipliers& operator=(Multipliers&& other) = default;

      void reset();
      [[nodiscard]] bool not_all_zero(size_t number_variables, double tolerance) const;
   };
} // namespace

#endif // UNO_MULTIPLIERS_H
