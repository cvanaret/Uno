#ifndef INTERIORPOINT_H
#define INTERIORPOINT_H

#include "LocalApproximation.hpp"
#include "MA57Solver.hpp"

/*! \class LocalApproximation
* \brief Local approximation
*
*  Local appromination of a nonlinear optimization problem (virtual class) 
*/
class InteriorPoint: public LocalApproximation {
	public:
		/*!
         *  Constructor
         */
		InteriorPoint();
		
		void allocate_solver(int number_variables, int number_constraints);
		
		LocalSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double objective_multiplier, double radius);
		
		LocalSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, const std::vector<double>& d, ConstraintPartition& constraint_partition, std::vector<double>& multipliers);
		
		LocalSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions);
		
		MA57Solver solver; /*!< Solver that solves the subproblem */
		
	private:
		std::vector<double> compute_slack_displacements(Problem& problem, std::vector<ConstraintType>& variable_status, int number_slacks, double mu, std::vector<double>& multipliers, std::vector<double>& s, std::vector<double>& solution);
};

#endif // INTERIORPOINT_H
