#ifndef LOCALSOLUTION_H
#define LOCALSOLUTION_H

#include <ostream>
#include <vector>
#include "Utils.hpp"
#include "Constraint.hpp"
#include "Phase.hpp"

/* see bqpd.f */
enum Status {OPTIMAL = 0,
			UNBOUNDED,
			BOUND_INCONSISTENCY,
			INFEASIBLE,
			INCORRECT_PARAMETER,
			LP_INSUFFICIENT_SPACE,
			HESSIAN_INSUFFICIENT_SPACE,
			SPARSE_INSUFFICIENT_SPACE,
			MAX_RESTARTS_REACHED,
			UNDEFINED};

/*! \struct LocalSolution
* \brief Solution of a local subproblem
*
*  Description of a local solution
*/
class LocalSolution {
	public:
		LocalSolution(std::vector<double>& x, int n, int m);
		
		Status status; /*!< Status of the computed step */
		Phase phase;
		std::vector<double> x; /*!< Primal variables in \f$\mathbf{R}^n\f$ */
		double norm; /*!< Norm of \f$x\f$ */
		double objective; /*!< Objective value */
		std::vector<double> multipliers; /*!< Approximate Lagrange multipliers/dual variables of the subproblem */
		ConstraintActivity active_set; /*!< Active set */
		ConstraintPartition constraint_partition; /*!< Partition of feasible and infeasible constraints */
		
		friend std::ostream& operator<< (std::ostream &stream, LocalSolution& step);
};

#endif // LOCALSOLUTION_H
