#ifndef TUBESTRATEGY_H
#define TUBESTRATEGY_H

#include "TwoPhaseStrategy.hpp"
#include "Tube.hpp"

/*! \class TubeStrategy
* \brief Step acceptance strategy based on a tube
*
*  Strategy that accepts or declines a trial step
*/
class TubeStrategy: public TwoPhaseStrategy {
	public:
		/*!
         *  Constructor that takes an optimization problem and a set of constants
         * 
         * \param problem: optimization problem
         * \param constants: set of constants
         */
		TubeStrategy(StepComputation& step_computation, StepConstants& constants, Tolerances& tolerances);
	
		Tube tube_restoration; /*!< Tube for restoration phase */
		Tube tube_optimality; /*!< Tube for optimality phase */
		
		bool check_step(Problem& problem, Iterate& current_iterate, Step& step, Phase& phase);
		
		void initialize(Iterate& current_iterate);
	
	protected:
		void phase_transition(Problem& problem, Iterate& current_iterate, Iterate& trial_iterate, Step& step, Phase& phase, ConstraintPartition& constraint_partition);
};

#endif // TUBESTRATEGY_H
