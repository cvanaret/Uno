#ifndef GLOBALIZATIONMECHANISM_H
#define GLOBALIZATIONMECHANISM_H

#include "Problem.hpp"
#include "GlobalizationStrategy.hpp"

/*! \class GlobalizationMechanism
* \brief Step control strategy
*
*  Strategy that promotes global convergence
*/
class GlobalizationMechanism {
	public:
		/*!
         *  Constructor
         * 
         * \param direction_computation: strategy to compute a descent direction
         * \param step_accept: strategy to accept or reject a step
         */
		GlobalizationMechanism(GlobalizationStrategy& globalization_strategy, int max_iterations);
		virtual ~GlobalizationMechanism();
		
		/*!
         *  Compute the next iterate from a given point
         * (purely virtual)
         * 
         * \param problem: optimization problem
         * \param current_iterate: current point and its evaluations
         */
		virtual Iterate compute_iterate(Problem& problem, Iterate& current_iterate) = 0;
		
		void initialize(Problem& problem, Iterate& current_iterate);
		
		/* references to allow polymorphism */
		GlobalizationStrategy& globalization_strategy;  /*!< Strategy to accept or reject a step */
		int max_iterations; /*!< Maximum number of iterations */
		int number_iterations; /*!< Current number of iterations */
};

#endif // GLOBALIZATIONMECHANISM_H
