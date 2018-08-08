#ifndef INTERIORPOINT_H
#define INTERIORPOINT_H

#include "Subproblem.hpp"
#include "MA57Solver.hpp"

/*! \class InteriorPoint
* \brief Interior Point Method
*
*  Implementation of an Interior Point Method
*/
class InteriorPoint: public Subproblem {
	public:
		/*!
         *  Constructor
         */
		InteriorPoint();
		
		void initialize(Problem& problem, Iterate& current_iterate, int number_variables, int number_constraints, double radius);
		
		LocalSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double objective_multiplier, double radius);
		
		LocalSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, const std::vector<double>& d, ConstraintPartition& constraint_partition, std::vector<double>& multipliers);
		
		LocalSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions);
		
		MA57Solver solver; /*!< Solver that solves the subproblem */
		
		double mu;
		
		/* slacks and multipliers */
		std::vector<double> multipliers;
		std::vector<double> slacks;
		std::vector<ConstraintType> variable_status;
		
		/* constants */
		double tau;
		double default_multiplier;
		double default_slack;

	private:
		COOMatrix generate_kkt_matrix(Problem& problem, Iterate& current_iterate, double objective_multiplier, std::vector<double> original_multipliers);
		std::vector<double> generate_rhs(Problem& problem, Iterate& current_iterate, double objective_multiplier, std::vector<double>& original_multipliers, std::vector<double> variable_lb, std::vector<double> variable_ub);
		std::vector<double> compute_slack_displacements(Problem& problem, double mu, std::vector<double>& solution);
};

#endif // INTERIORPOINT_H
