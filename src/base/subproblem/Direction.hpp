#ifndef DIRECTION_H
#define DIRECTION_H

#include <ostream>
#include <vector>
#include <set>
#include <functional>
#include "Vector.hpp"
#include "Problem.hpp"
#include "Iterate.hpp"
#include "Constraint.hpp"
#include "Phase.hpp"

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
    Direction(std::vector<double>& x, Multipliers& multipliers);
    std::vector<double> x; /*!< Primal variables */
    Multipliers multipliers; /*!< Multipliers */
    double objective_multiplier; /*!< Objective multiplier */

    Status status; /*!< Status of the solution */
    Phase phase; /*!< Current phase */

    double norm; /*!< Norm of \f$x\f$ */
    double objective; /*!< Objective value */
    ActiveSet active_set; /*!< Active set */
    std::set<int> inactive_set; /*!< Inactive set */
    ConstraintPartition constraint_partition; /*!< Partition of feasible and infeasible constraints */
    
    // this function computes the predicted reduction of the direction for a given step length
    std::function<double(double step_length)> predicted_reduction;

    friend std::ostream& operator<<(std::ostream &stream, Direction& step);
};

#endif // DIRECTION_H
