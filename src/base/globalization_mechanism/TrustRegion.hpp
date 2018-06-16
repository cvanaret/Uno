#ifndef TRUSTREGION_H
#define TRUSTREGION_H

#include "GlobalizationMechanism.hpp"

/*! \class TrustRegion
* \brief Trust region
*
*  Trust region strategy
*/
class TrustRegion: public GlobalizationMechanism {
	public:
		/*!
         *  Constructor
         * 
         * \param direction_computation: strategy to compute a descent direction
         * \param step_accept: strategy to accept or reject a step
         * \param initial_radius: initial trust region radius
         */
		TrustRegion(GlobalizationStrategy& globalization_strategy, double initial_radius, int max_iterations = 100);
		
		/*!
         *  Compute the next iterate from a given point
         * 
         * \param problem: optimization problem
         * \param current_iterate: current point and its evaluations
         */
		Iterate compute_iterate(Problem& problem, Iterate& current_iterate);
		
		double radius; /*!< Current trust region radius */

	private:
		void correct_multipliers(Problem& problem, LocalSolution& solution);
		/*!
         *  Determine if a new iterate has been found
         * 
         * \param success: true if a new step has been accepted
         * \param iteration: current iteration number
         * \param radius: current trust region radius
         */
		bool termination(bool success, int iteration, double radius);
	
		/*!
         *  Keep track of min/max/average radius
         * 
         * \param radius: current trust region radius
         */
		void record_radius(double radius);

		double radius_max_; /*!< Maximum used radius */
		double radius_min_; /*!< Minimum used radius */
		double radius_sum_; /*!< Sum of all used radii */
		
		double activity_tolerance_;
};

#endif // TRUSTREGION_H
