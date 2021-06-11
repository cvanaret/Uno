#ifndef GLOBALIZATIONMECHANISM_H
#define GLOBALIZATIONMECHANISM_H

#include "Problem.hpp"
#include "GlobalizationStrategy.hpp"
#include "Statistics.hpp"

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
   GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations);
   virtual ~GlobalizationMechanism() = default;

   virtual Iterate initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) = 0;
   virtual std::pair<Iterate, Direction> compute_acceptable_iterate(Statistics& statistics, const Problem& problem, Iterate& current_iterate) = 0;

   /* references to allow polymorphism */
   ConstraintRelaxationStrategy& relaxation_strategy;
   int max_iterations; /*!< Maximum number of iterations */
   int number_iterations; /*!< Current number of iterations */

protected:
   // preallocated vector to receive the trial primal variables
   std::vector<double> trial_primals_;

   Iterate assemble_trial_iterate(const Problem& problem, const Iterate& current_iterate, const Direction& direction, double step_length);
   virtual void print_acceptance_(const Iterate& iterate);
   static void print_warning_(const char* message);
};

#endif // GLOBALIZATIONMECHANISM_H
