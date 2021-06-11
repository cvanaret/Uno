#ifndef GLOBALIZATIONSTRATEGY_H
#define GLOBALIZATIONSTRATEGY_H

#include <cmath>
#include <optional>
#include "ConstraintRelaxationStrategy.hpp"
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
   GlobalizationStrategy() = default;
   virtual ~GlobalizationStrategy() = default;

   virtual void initialize(Statistics& statistics, const Iterate& first_iterate) = 0;
   virtual bool check_acceptance(Statistics& statistics, ProgressMeasures& current_progress, ProgressMeasures& trial_progress,
         double objective_multiplier, double predicted_reduction) = 0;
   virtual void reset() = 0;
   virtual void notify(Iterate& current_iterate) = 0;
};

#endif // GLOBALIZATIONSTRATEGY_H
