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
		
		void initialize(Problem& problem, Iterate& current_iterate, int number_variables, int number_constraints, bool use_trust_region);
		
		LocalSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double objective_multiplier, double radius);
		
		LocalSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, const std::vector<double>& d, ConstraintPartition& constraint_partition, std::vector<double>& multipliers);
		
		LocalSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions);
		
		MA57Solver solver; /*!< Solver that solves the subproblem */
		double tau;
		double default_multiplier;
		double default_slack;
		
	private:
		std::vector<double> compute_slack_displacements(Problem& problem, std::vector<ConstraintType>& variable_status, double mu, std::vector<double>& multipliers, std::vector<double>& slacks, std::vector<double>& solution);
};

#endif // INTERIORPOINT_H
