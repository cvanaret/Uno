// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DIRECTION_H
#define UNO_DIRECTION_H

#include <vector>
#include <optional>
#include <ostream>
#include "ingredients/subproblems/SubproblemStatus.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Multipliers.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   struct ActiveConstraints {
      std::vector<size_t> at_lower_bound; /*!< List of constraint indices at their lower bound */
      std::vector<size_t> at_upper_bound; /*!< List of constraint indices at their upper bound */

      explicit ActiveConstraints(size_t capacity) {
         this->at_lower_bound.reserve(capacity);
         this->at_upper_bound.reserve(capacity);
      }
   };

   class Direction {
   public:
      Direction(size_t number_variables, size_t number_constraints);

      size_t number_variables;
      size_t number_constraints;

      Vector<double> primals; /*!< Primal variables */
      Multipliers multipliers; /*!< Multipliers */
      Multipliers feasibility_multipliers; /*!< Multipliers */

      SubproblemStatus status{SubproblemStatus::OPTIMAL}; /*!< Status of the solution */

      double norm{INF<double>}; /*!< Norm of \f$x\f$ */
      double subproblem_objective{INF<double>}; /*!< Objective value */
      ActiveConstraints active_bounds; /*!< Active bound constraints */

      void set_dimensions(size_t new_number_variables, size_t new_number_constraints);
      void reset();

      friend std::ostream& operator<<(std::ostream& stream, const Direction& direction);
   };
} // namespace

#endif // UNO_DIRECTION_H
