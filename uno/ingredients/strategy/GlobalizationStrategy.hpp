#ifndef UNO_GLOBALIZATIONSTRATEGY_H
#define UNO_GLOBALIZATIONSTRATEGY_H

#include "optimization/Iterate.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "tools/Statistics.hpp"

/*! \class GlobalizationStrategy
 * \brief Step acceptance strategy
 *
 *  Strategy that accepts or declines a trial step (virtual class)
 */
class GlobalizationStrategy {
public:
   GlobalizationStrategy() = default;
   virtual ~GlobalizationStrategy() = default;

   virtual void initialize(Statistics& statistics, const Iterate& first_iterate) = 0;
   virtual bool check_acceptance(Statistics& statistics, const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress,
         double objective_multiplier, double predicted_reduction) = 0;
   virtual void reset() = 0;
   virtual void notify(Iterate& current_iterate) = 0;
};

#endif // UNO_GLOBALIZATIONSTRATEGY_H
