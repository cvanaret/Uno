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
    //SubproblemSolution(std::vector<double>& x, Multipliers& multipliers, Status& status, Phase& phase, bool phase_1_required, double norm, double objective, ActiveSet& active_set, ConstraintPartition& constraint_partition);
    SubproblemSolution(std::vector<double>& x, Multipliers& multipliers, ActiveSet& active_set, ConstraintPartition& constraint_partition);
    std::vector<double> x; /*!< Primal variables */
    Multipliers multipliers; /*!< Multipliers */
    
    Status status; /*!< Status of the solution */
    Phase phase; /*!< Current phase */
    bool phase_1_required;

    double norm; /*!< Norm of \f$x\f$ */
    double objective; /*!< Objective value */
    //ObjectiveTerms objective_terms; /*!< Decomposition of the objective value in quadratic and linear terms */
    ActiveSet active_set; /*!< Active set */
    ConstraintPartition constraint_partition; /*!< Partition of feasible and infeasible constraints */

    friend std::ostream& operator<<(std::ostream &stream, SubproblemSolution& step);
};

#endif // SUBPROBLEMSOLUTION_H
