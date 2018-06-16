#ifndef TRUSTLINESEARCH_H
#define TRUSTLINESEARCH_H

#include "GlobalizationMechanism.hpp"

/*! \class LineSearch
* \brief Line-search
*
*  Line-search strategy
*/
class TrustLineSearch: public GlobalizationMechanism {
	public:
		/*!
         *  Constructor
         */
		TrustLineSearch(GlobalizationStrategy& globalization_strategy, double initial_radius, int max_iterations = 30, double ratio = 0.5);
		
		/*!
         *  Compute the next iterate from a given point
         * 
         * \param problem: optimization problem
         * \param current_iterate: current point and its evaluations
         */
		Iterate compute_iterate(Problem& problem, Iterate& current_iterate);
		
		/* ratio of step length update in ]0, 1[ */
		double ratio;
		double radius; /*!< Current trust region radius */
		
	private:
		bool termination(bool is_accepted, int iteration);
		
		void correct_multipliers(Problem& problem, LocalSolution& solution);
		
		double activity_tolerance_;
};

#endif // TRUSTLINESEARCH_H
