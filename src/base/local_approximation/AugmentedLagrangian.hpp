#ifndef AUGMENTEDLAGRANGIAN_H
#define AUGMENTEDLAGRANGIAN_H

#include "StepComputation.hpp"
#include "QPSolver.hpp"

/*! \class AugmentedLagrangian
* \brief Augmented Lagrangian strategy
*
*  Strategy that computes a descent direction based on an augmented Lagrangian
*/
class AugmentedLagrangian: public StepComputation {
	public:
		/*!
         *  Constructor
         * 
         * \param solver: solver that computes the step
         */
		AugmentedLagrangian(QPSolver& solver);
		
		/*!
         *  Compute the descent direction for a given point and phase, with an zero initial solution
         * 
         * \param problem: optimization problem
         * \param phase: current phase (optimality or feasibility restoration)
         * \param current_point: current point and its evaluations
         * \param radius: possible radius of the trust region
         */
		Step compute_optimality_direction(Problem& problem, Iterate& current_point, double radius);
		
		/*!
         *  Compute the descent direction for a given point and phase, with a given initial solution
         * 
         * \param problem: optimization problem
         * \param phase: current phase (optimality or feasibility restoration)
         * \param current_point: current point and its evaluations
         * \param radius: possible radius of the trust region
         * \param d0: initial solution
         */
		Step compute_infeasibility_direction(Problem& problem, ConstraintPartition& constraint_partition, Iterate& current_point, double radius, std::vector<double>& d0);
		
		/* use a reference to allow polymorphism */
		QPSolver& solver; /*!< Solver that computes the step */
};

#endif // AUGMENTEDLAGRANGIAN_H
