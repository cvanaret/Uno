#ifndef LPAPPROXIMATION_H
#define LPAPPROXIMATION_H

#include "LocalApproximation.hpp"
#include "LPSolver.hpp"

/*! \class LPApproximation
* \brief LP local approximation
*
*  Linear approximation
*/
class LPApproximation: public LocalApproximation {
	public:
		/*!
         *  Constructor
         * 
         * \param solver: solver that solves the subproblem
         */
		LPApproximation(LPSolver& solver);
		
		/*!
         *  Compute the descent direction for a given point and phase
         * 
         * \param problem: optimization problem
         * \param phase: current phase (optimality or feasibility restoration)
         * \param current_iterate: current point and its evaluations
         */
		LocalSolution compute_direction(Problem& problem, Iterate& current_iterate, Phase& phase);
		
		// TODO static list of available solvers (BQPD)
		
		/* use a reference to allow polymorphism */
		LPSolver& solver; /*!< Solver that solves the subproblem */
	
	private:
		/*!
         *  Generate a linear local approximation
         * 
         * \param problem: optimization problem
         * \param current_iterate: current point and its evaluations
         */
		LP generate_LP(Problem& problem, Iterate& current_iterate) const;
};

#endif // LPAPPROXIMATION_H
