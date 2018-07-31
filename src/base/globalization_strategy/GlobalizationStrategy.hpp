#ifndef GLOBALIZATIONSTRATEGY_H
#define GLOBALIZATIONSTRATEGY_H

#include "Problem.hpp"
#include "LocalApproximation.hpp"
#include "Iterate.hpp"
#include "LocalSolution.hpp"
#include "Constraint.hpp"

/*! \class GlobalizationStrategy
* \brief Step acceptance strategy
*
*  Strategy that accepts or declines a trial step (virtual class)
*/
class GlobalizationStrategy {
	public:
		/*!
         *  Constructor that takes an optimization problem and a tolerance
         * 
         * \param problem: optimization problem
         * \param constants: set of constants
         */
		GlobalizationStrategy(LocalApproximation& local_approximation, double tolerance);
		virtual ~GlobalizationStrategy();
		
		LocalApproximation& local_approximation;
		double tolerance; /*!< Tolerance of the termination criteria */
		
		virtual LocalSolution compute_step(Problem& problem, Iterate& current_iterate, double radius) = 0;
		
		/*!
         *  Check the validity of a step
         *  Purely virtual method (only implemented in subclasses)
         */
		virtual bool check_step(Problem& problem, Iterate& current_iterate, LocalSolution& solution, double step_length = 1.) = 0;
		
		virtual void initialize(Problem& problem, Iterate& current_iterate) = 0;
		
		std::vector<double> compute_lagrangian_gradient(Problem& problem, Iterate& current_iterate, double objective_multiplier, std::vector<double>& multipliers);
		
		virtual double compute_KKT_error(Problem& problem, Iterate& current_iterate) = 0;
		
		double compute_complementarity_error(const Problem& problem, Iterate& current_iterate);
};

#endif // GLOBALIZATIONSTRATEGY_H
