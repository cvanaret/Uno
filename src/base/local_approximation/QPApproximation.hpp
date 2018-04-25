#ifndef QPAPPROXIMATION_H
#define QPAPPROXIMATION_H

#include "LocalApproximation.hpp"
#include "QPSolver.hpp"

/*! \class QPApproximation
* \brief QP local approximation
*
*  Quadratic approximation
*/
class QPApproximation: public LocalApproximation {
	public:
		/*!
         *  Constructor
         * 
         * \param solver: solver that solves the subproblem
         */
		QPApproximation(QPSolver& solver);
		
		void allocate_solver(int number_variables, int number_constraints);
		
		LocalSolution compute_optimality_step(Problem& problem, Iterate& current_point, double objective_multiplier, double radius);
		
		LocalSolution compute_infeasibility_step(Problem& problem, Iterate& current_point, double radius, const std::vector<double>& d, ConstraintPartition& constraint_partition, std::vector<double>& multipliers);
		
		LocalSolution compute_l1_penalty_step(Problem& problem, Iterate& current_point, double radius, double penalty_parameter, PenaltyConstraints penalty_constraints);
		
		/* use a reference to allow polymorphism */
		QPSolver& solver; /*!< Solver that solves the subproblem */

	private:
		/*!
         *  Generate a quadratic local approximation
         * 
         * \param problem: optimization problem
         * \param current_point: current point and its evaluations
         */
		QP generate_optimality_qp(Problem& problem, Iterate& current_point, double objective_multiplier, double radius);
		
		QP generate_infeasibility_qp(Problem& problem, Iterate& current_point, double radius, ConstraintPartition& constraint_partition, std::vector<double>& multipliers);
		
		QP generate_qp(Problem& problem, Iterate& current_point, double radius);
		
		QP generate_l1_penalty_qp(Problem& problem, Iterate& current_point, double radius, double penalty_parameter, PenaltyConstraints penalty_constraints);
		
		void set_constraints(Problem& problem, QP& qp, Iterate& current_point);
		
		void set_optimality_objective(Problem& problem, QP& qp, Iterate& current_point);
		
		void set_infeasibility_objective(Problem& problem, QP& qp, Iterate& current_point, ConstraintPartition& constraint_partition);
		
		/*!
         *  Generate a quadratic local approximation
         * 
         * \param problem: optimization problem
         * \param current_point: current point and its evaluations
         */
		//QP generate_EQP(Problem& problem, Iterate& current_point);
		
		/*!
         *  Generate a quadratic local approximation
         * 
         * \param problem: optimization problem
         * \param current_point: current point and its evaluations
         */
		//QP generate_EQPI(Problem& problem, Iterate& current_point);
};

#endif // QPAPPROXIMATION_H
