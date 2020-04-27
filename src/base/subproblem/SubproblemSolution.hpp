#ifndef SUBPROBLEMSOLUTION_H
#define SUBPROBLEMSOLUTION_H

#include <ostream>
#include <vector>
#include "Utils.hpp"
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

struct ObjectiveTerms {
    double linear;
    double quadratic;
};

/*! \struct SubproblemSolution
 * \brief Solution of a local subproblem
 *
 *  Description of a local solution
 */
class SubproblemSolution {
public:
    SubproblemSolution(std::vector<double>& x, Multipliers& multipliers);
    std::vector<double> x; /*!< Primal variables */
    Multipliers multipliers; /*!< Multipliers */
    double objective_multiplier; /*!< Objective multiplier */
    
    Status status; /*!< Status of the solution */
    Phase phase; /*!< Current phase */
    bool phase_1_required;

    double norm; /*!< Norm of \f$x\f$ */
    double objective; /*!< Objective value */
    bool is_descent_direction;
    ActiveSet active_set; /*!< Active set */
    ConstraintPartition constraint_partition; /*!< Partition of feasible and infeasible constraints */

    friend std::ostream& operator<<(std::ostream &stream, SubproblemSolution& step);
};

#endif // SUBPROBLEMSOLUTION_H
