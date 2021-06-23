#ifndef DIRECTION_H
#define DIRECTION_H

#include <ostream>
#include <vector>
#include <set>
#include <functional>
#include <cmath>
#include "Vector.hpp"
#include "Problem.hpp"
#include "Iterate.hpp"
#include "Constraint.hpp"

/* see bqpd.f */
enum Status {
   OPTIMAL = 0,
   UNBOUNDED_PROBLEM,
   BOUND_INCONSISTENCY,
   INFEASIBLE,
   INCORRECT_PARAMETER,
   LP_INSUFFICIENT_SPACE,
   HESSIAN_INSUFFICIENT_SPACE,
   SPARSE_INSUFFICIENT_SPACE,
   MAX_RESTARTS_REACHED,
   UNDEFINED
};

/*! \struct SubproblemSolution
 * \brief Solution of a local subproblem
 *
 *  Description of a local solution
 */
class Direction {
public:
   Direction(size_t number_variables, size_t number_constraints);
   Direction(std::vector<double>& x, Multipliers& multipliers);
   std::vector<double> x; /*!< Primal variables */
   Multipliers multipliers; /*!< Multipliers */
   double objective_multiplier{1.}; /*!< Objective multiplier */

   Status status{OPTIMAL}; /*!< Status of the solution */
   bool is_relaxed{false};

   double norm{0.}; /*!< Norm of \f$x\f$ */
   double objective{INFINITY}; /*!< Objective value */
   ActiveSet active_set; /*!< Active set */
   std::vector<int> inactive_set; /*!< Inactive set */
   ConstraintPartition constraint_partition; /*!< Partition of feasible and infeasible constraints */

   // this function computes the predicted reduction of the direction for a given step length
   std::function<double(double step_length)> predicted_reduction{nullptr};

   friend std::ostream& operator<<(std::ostream& stream, const Direction& step);
};

#endif // DIRECTION_H
