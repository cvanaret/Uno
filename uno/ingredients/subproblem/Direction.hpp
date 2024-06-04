// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DIRECTION_H
#define UNO_DIRECTION_H

#include <vector>
#include <optional>
#include <ostream>
#include "SubproblemStatus.hpp"
#include "optimization/Multipliers.hpp"
#include "tools/Infinity.hpp"

// forward declaration
template <typename ElementType>
class Vector;

/*! \struct ConstraintActivity
* \brief Constraints at lower or upper bound at the optimum solution
*
*  Description of the active or infeasible constraints: at lower or upper bound at the optimum solution
*/
struct ActiveConstraints {
   std::vector<size_t> at_lower_bound; /*!< List of constraint indices at their lower bound */
   std::vector<size_t> at_upper_bound; /*!< List of constraint indices at their upper bound */

   explicit ActiveConstraints(size_t capacity) {
      this->at_lower_bound.reserve(capacity);
      this->at_upper_bound.reserve(capacity);
   }
};

struct ActiveSet {
   ActiveConstraints constraints; /*!< List of general constraints */
   ActiveConstraints bounds; /*!< List of bound constraints */

   ActiveSet(size_t number_variables, size_t number_constraints): constraints(number_constraints), bounds(number_variables) { }
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
   ActiveSet active_set; /*!< Active set */

   void set_dimensions(size_t new_number_variables, size_t new_number_constraints);
   void reset();

   friend std::ostream& operator<<(std::ostream& stream, const Direction& direction);
};

#endif // UNO_DIRECTION_H
