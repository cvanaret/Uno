#ifndef GLOBALIZATIONSTRATEGY_H
#define GLOBALIZATIONSTRATEGY_H

#include <cmath>
#include <optional>
#include "FeasibilityStrategy.hpp"
#include "Problem.hpp"
#include "Subproblem.hpp"
#include "Iterate.hpp"
#include "Direction.hpp"
#include "Statistics.hpp"

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
   explicit GlobalizationStrategy(FeasibilityStrategy& feasibility_strategy, Subproblem& subproblem);
   virtual ~GlobalizationStrategy() = default;

   FeasibilityStrategy& feasibility_strategy;
   Subproblem& subproblem;

   virtual Iterate initialize(Statistics& statistics, Problem& problem, std::vector<double>& x, Multipliers& multipliers) = 0;
   virtual std::optional<Iterate> check_acceptance(Statistics& statistics, Problem& problem, Iterate& current_iterate, Direction& direction,
         double step_length) = 0;

protected:
   // preallocated vector to receive the trial primal variables
   std::vector<double> trial_primals_;
};

#endif // GLOBALIZATIONSTRATEGY_H
