#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <vector>

enum ConstraintType { EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED };

/*! \struct ConstraintActivity
* \brief Constraints at lower or upper bound at the optimum solution
*
*  Description of the active or infeasible constraints: at lower or upper bound at the optimum solution
*/
struct ActiveConstraints {
   std::vector<int> at_lower_bound; /*!< List of constraint indices at their lower bound */
   std::vector<int> at_upper_bound; /*!< List of constraint indices at their upper bound */
};

struct ActiveSet {
   ActiveConstraints constraints; /*!< List of general constraints */
   ActiveConstraints bounds; /*!< List of bound constraints */
};

enum ConstraintFeasibility { FEASIBLE, INFEASIBLE_LOWER, INFEASIBLE_UPPER };

struct ConstraintPartition {
   std::vector<int> feasible; /*!< Indices of the feasible constraints */
   std::vector<int> infeasible; /*!< Indices of the infeasible constraints */
   std::vector<ConstraintFeasibility> constraint_feasibility;

   explicit ConstraintPartition(size_t number_constraints);
};

struct Range {
   double lb;
   double ub;
};

struct Multipliers {
   std::vector<double> lower_bounds; /*!< Multipliers of the lower bound constraints */
   std::vector<double> upper_bounds; /*!< Multipliers of the lower bound constraints */
   std::vector<double> constraints; /*!< Multipliers of the general constraints */
   double objective{1.};

   Multipliers() = default;
   Multipliers(size_t number_variables, size_t number_constraints);
};

#endif // CONSTRAINT_H
