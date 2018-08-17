#ifndef QPAPPROXIMATION_H
#define QPAPPROXIMATION_H

#include "Subproblem.hpp"
#include "QPSolver.hpp"

/*! \class QPApproximation
* \brief QP local approximation
*
*  Quadratic approximation
*/
class QPApproximation: public Subproblem {
	public:
		/*!
         *  Constructor
         * 
         * \param solver: solver that solves the subproblem
         */
		QPApproximation(QPSolver& solver);
		
		void initialize(Problem& problem, Iterate& current_iterate, int number_variables, int number_constraints, bool use_trust_region);
		
		LocalSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double radius);
		
		LocalSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, const std::vector<double>& d, ConstraintPartition& constraint_partition, std::vector<double>& multipliers);
		
		LocalSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions);
		
		/* use a reference to allow polymorphism */
		QPSolver& solver; /*!< Solver that solves the subproblem */

	private:
		/*!
         *  Generate a quadratic local approximation
         * 
         * \param problem: optimization problem
         * \param current_iterate: current point and its evaluations
         */
		QP generate_optimality_qp_(Problem& problem, Iterate& current_iterate, double radius);
		
		QP generate_infeasibility_qp_(Problem& problem, Iterate& current_iterate, double radius, ConstraintPartition& constraint_partition, std::vector<double>& multipliers);
		
		QP generate_qp_(Problem& problem, Iterate& current_iterate, double radius);
		
		QP generate_l1_penalty_qp_(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions);
		
		void set_constraints_(Problem& problem, QP& qp, Iterate& current_iterate);
		
		void set_optimality_objective_(Problem& problem, QP& qp, Iterate& current_iterate);
		
		void set_infeasibility_objective_(Problem& problem, QP& qp, Iterate& current_iterate, ConstraintPartition& constraint_partition);
		
		/*!
         *  Generate a quadratic local approximation
         * 
         * \param problem: optimization problem
         * \param current_iterate: current point and its evaluations
         */
		//QP generate_EQP(Problem& problem, Iterate& current_iterate);
		
		/*!
         *  Generate a quadratic local approximation
         * 
         * \param problem: optimization problem
         * \param current_iterate: current point and its evaluations
         */
		//QP generate_EQPI(Problem& problem, Iterate& current_iterate);
};

#endif // QPAPPROXIMATION_H
