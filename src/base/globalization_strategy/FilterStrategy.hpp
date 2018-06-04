#ifndef FILTERSTRATEGY_H
#define FILTERSTRATEGY_H

#include "TwoPhaseStrategy.hpp"
#include "Filter.hpp"

/*! \class FilterStrategy
* \brief Step acceptance strategy based on a filter
*
*  Strategy that accepts or declines a trial step
*/
class FilterStrategy: public TwoPhaseStrategy {
	public:
		/*!
         *  Constructor that takes an optimization problem, filters for restoration and optimality, and a set of constants
         */
		FilterStrategy(LocalApproximation& local_approximation, Filter& filter_restoration, Filter& filter_optimality, LocalSolutionConstants& constants, Tolerances& tolerances, double tolerance);
	
		/* use references to allow polymorphism */
		Filter& filter_restoration; /*!< Filter for the restoration phase */
		Filter& filter_optimality; /*!< Filter for the optimality phase */
		Tolerances& tolerances; /*!< Tolerances */
		
		/*!
         *  Check the validity of a step
         *  Implements the purely virtual method of the superclass
         */
		bool check_step(Problem& problem, Iterate& current_iterate, LocalSolution& solution, double step_length = 1.);
		
		void initialize(Problem& problem, Iterate& current_iterate);
		
	private:
		OptimalityStatus compute_status(Problem& problem, Iterate& current_iterate, double step_norm);
};

#endif // FILTERSTRATEGY_H
