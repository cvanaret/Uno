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
		
		void initialize(Problem& problem, Iterate& current_iterate, int number_variables, int number_constraints, bool use_trust_region);
		
		LocalSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double radius);
		
		LocalSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, const std::vector<double>& d, ConstraintPartition& constraint_partition, std::vector<double>& multipliers);
		
		LocalSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions);
		
		MA57Solver solver; /*!< Solver that solves the subproblem */
		
		double mu;
		std::vector<ConstraintType> variable_status;
		
		/* slacks and multipliers */
		std::vector<double> bound_multipliers;
		std::vector<double> constraint_multipliers;
		int number_slacks;
		
		/* constants */
		double tau_min;
		double default_multiplier;
		double k_sigma;
		double inertia_term = 0.;
		
		MA57Data data;

	private:
		double project_variable_in_bounds(double current_value, double lb, double ub);
		std::vector<double> estimate_initial_multipliers(Problem& problem, Iterate& current_iterate);
		COOMatrix generate_kkt_matrix(Problem& problem, Iterate& current_iterate, std::vector<double>& original_multipliers, std::vector<double>& variable_lb, std::vector<double>& variable_ub);
		std::vector<double> generate_rhs(Problem& problem, Iterate& current_iterate, std::vector<double>& original_multipliers, std::vector<double>& variable_lb, std::vector<double>& variable_ub);
		std::vector<double> compute_slack_displacements(Problem& problem, Iterate& current_iterate, std::vector<double>& solution);
		std::vector<double> compute_bound_multiplier_displacements(Problem& problem, Iterate& current_iterate, std::vector<double>& solution, std::vector<double>& variable_lb, std::vector<double>& variable_ub);
		double update_barrier_parameter(Problem& problem, Iterate& current_iterate);
};

#endif // INTERIORPOINT_H
