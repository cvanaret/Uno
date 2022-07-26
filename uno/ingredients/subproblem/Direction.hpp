// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DIRECTION_H
#define UNO_DIRECTION_H

#include <vector>
#include <optional>
#include <ostream>
#include "optimization/Multipliers.hpp"
#include "tools/Infinity.hpp"

enum class Status {
   OPTIMAL = 0,
   UNBOUNDED_PROBLEM,
   INFEASIBLE
};

/*! \struct ConstraintActivity
* \brief Constraints at lower or upper bound at the optimum solution
*
*  Description of the active or infeasible constraints: at lower or upper bound at the optimum solution
*/
struct ActiveConstraints {
   std::vector<size_t> at_lower_bound{}; /*!< List of constraint indices at their lower bound */
   std::vector<size_t> at_upper_bound{}; /*!< List of constraint indices at their upper bound */
};

struct ActiveSet {
   ActiveConstraints constraints{}; /*!< List of general constraints */
   ActiveConstraints bounds{}; /*!< List of bound constraints */
};

struct ConstraintPartition {
   std::vector<size_t> feasible{}; /*!< Indices of the feasible constraints */
   std::vector<size_t> infeasible{}; /*!< Indices of the infeasible constraints */
   std::vector<size_t> lower_bound_infeasible{}; /*!< Indices of the lower-bound infeasible constraints */
   std::vector<size_t> upper_bound_infeasible{}; /*!< Indices of the upper_bound infeasible constraints */

   explicit ConstraintPartition(size_t number_constraints);
};

class Direction {
public:
   Direction(size_t number_variables, size_t number_constraints);

   std::vector<double> primals; /*!< Primal variables */
   Multipliers multipliers; /*!< Multipliers */
   double objective_multiplier{1.}; /*!< Objective multiplier */

   Status status{Status::OPTIMAL}; /*!< Status of the solution */

   double norm{INF}; /*!< Norm of \f$x\f$ */
   double objective{INF}; /*!< Objective value */
   ActiveSet active_set{}; /*!< Active set */
   std::optional<ConstraintPartition> constraint_partition{std::nullopt}; /*!< Optional partition of feasible and infeasible constraints */

   friend std::ostream& operator<<(std::ostream& stream, const Direction& step);
};

#endif // UNO_DIRECTION_H
