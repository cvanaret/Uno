#ifndef TUBE_H
#define TUBE_H

#include "Step.hpp"

/*! \class Tube
* \brief Tube
*
*  Definition of a tube
*/
class Tube {
	public:
		/*!
         *  Constructor that takes an upper bound and a set of constants
         * 
         * \param upper_bound: upper bound of the tube
         * \param constants: set of constants
         */
		Tube(double upper_bound, StepConstants& constants);
		
		double upper_bound; /*!< Upper bound */
		StepConstants constants; /*!< Set of constants */
		
		/*!
         *  Update the upper bound
         * 
         * \param current_constraint_residual: current residual of the constraints
         * \param new_constraint_residual: residual of the constraints of the trial point
         */
		void update(double current_constraint_residual, double trial_constraint_residual);
		
		/*!
         *  Check whether new constraint violation is acceptable
         * 
         * \param constraint_residual: residual of the constraints
         */
		bool query(double constraint_residual);

	protected:

};

#endif // TUBE_H
