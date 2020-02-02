#ifndef LOCALSOLUTION_H
#define LOCALSOLUTION_H

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

/*! \struct LocalSolution
 * \brief Solution of a local subproblem
 *
 *  Description of a local solution
 */
class LocalSolution {
    public:
        LocalSolution(std::vector<double>& x, std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers);

        Status status; /*!< Status of the solution */
        Phase phase; /*!< Current phase */
        bool phase_1_required;
        std::vector<double> x; /*!< Primal variables */
        std::vector<double> bound_multipliers; /*!< Multipliers of the bound constraints */
        std::vector<double> constraint_multipliers; /*!< Multipliers of the general constraints */
        double norm; /*!< Norm of \f$x\f$ */
        double objective; /*!< Objective value */
        ObjectiveTerms objective_terms; /*!< Decomposition of the objective value in quadratic and linear terms */
        ConstraintActivity active_set; /*!< Active set */
        ConstraintPartition constraint_partition; /*!< Partition of feasible and infeasible constraints */

        friend std::ostream& operator<<(std::ostream &stream, LocalSolution& step);
};

#endif // LOCALSOLUTION_H
