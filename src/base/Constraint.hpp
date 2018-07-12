#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <vector>

enum ConstraintType {EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES};

struct PenaltyDimensions {
	int number_additional_variables;
	int number_constraints;
};

/*! \struct ConstraintActivity
* \brief Constraints at lower or upper bound at the optimum solution
*
*  Description of the active or infeasible constraints: at lower or upper bound at the optimum solution
*/
struct ConstraintActivity {
	std::vector<int> at_lower_bound; /*!< List of constraint indices at their lower bound */
	std::vector<int> at_upper_bound; /*!< List of constraint indices at their upper bound */
};

enum ConstraintFeasibility {FEASIBLE, INFEASIBLE_LOWER, INFEASIBLE_UPPER};

struct ConstraintPartition {
	std::vector<int> feasible_set; /*!< Indices of the feasible constraints */
	std::vector<int> infeasible_set; /*!< Indices of the infeasible constraints */
	std::vector<ConstraintFeasibility> status;
};

#endif // CONSTRAINT_H
