#ifndef LOCALAPPROXIMATION_H
#define LOCALAPPROXIMATION_H

#include <vector>
#include "Problem.hpp"
#include "Iterate.hpp"
#include "Phase.hpp"
#include "LocalSolution.hpp"
#include "Constraint.hpp"

/*! \class LocalApproximation
* \brief Local approximation
*
*  Local appromination of a nonlinear optimization problem (virtual class) 
*/
class LocalApproximation {
	public:
		/*!
         *  Constructor
         * 
         * \param solver: solver that solves the subproblem
         * \param name: name of the strategy
         */
		LocalApproximation(std::string name);
		
		virtual LocalSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double objective_multiplier, double radius) = 0;
		
		virtual LocalSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, const std::vector<double>& d, ConstraintPartition& constraint_partition, std::vector<double>& multipliers) = 0;
		
		virtual LocalSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyConstraints penalty_constraints) = 0;

		virtual void allocate_solver(int number_variables, int number_constraints) = 0;

		std::string name; /*!< Name of the strategy */
		int number_subproblems_solved;
};

#endif // LOCALAPPROXIMATION_H
