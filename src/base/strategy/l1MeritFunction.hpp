#ifndef PENALTYMERITFUNCTION_H
#define PENALTYMERITFUNCTION_H

#include <vector>
#include "GlobalizationStrategy.hpp"
#include "ConstraintRelaxationStrategy.hpp"
#include "Constraint.hpp"

class l1MeritFunction : public GlobalizationStrategy {
public:
   /*!
    *  Constructor that takes an optimization problem and a set of constants
    */
   explicit l1MeritFunction();

   void initialize(Statistics& statistics, const Iterate& first_iterate) override;
   bool check_acceptance(Statistics& statistics, ProgressMeasures& current_progress, ProgressMeasures& trial_progress, double objective_multiplier,
         double predicted_reduction) override;
   void reset() override;
   void notify(Iterate& current_iterate) override;

private:
   double decrease_fraction_;
};

#endif // PENALTYMERITFUNCTION_H
