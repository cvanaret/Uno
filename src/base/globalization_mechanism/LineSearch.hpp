#ifndef LINESEARCH_H
#define LINESEARCH_H

#include "GlobalizationMechanism.hpp"

/*! \class LineSearch
* \brief Line-search
*
*  Line-search strategy
*/
class LineSearch: public GlobalizationMechanism {
	public:
		/*!
         *  Constructor
         */
		LineSearch(GlobalizationStrategy& globalization_strategy, int max_iterations = 30, double ratio = 0.5);
		
		/*!
         *  Compute the next iterate from a given point
         * 
         * \param problem: optimization problem
         * \param current_iterate: current point and its evaluations
         */
		Iterate compute_iterate(Problem& problem, Iterate& current_iterate);
		
		double step_length;
		/* ratio of step length update in ]0, 1[ */
		double ratio;
		
	private:
		bool termination(bool is_accepted, int iteration);
};

#endif // LINESEARCH_H
