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
   explicit l1MeritFunction(Subproblem& subproblem);

   void initialize(Statistics& statistics, const Iterate& first_iterate) override;
   bool check_acceptance(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, Direction& direction, double step_length) override;
   void reset() override;
   void notify(Iterate& current_iterate) override;

private:
   double decrease_fraction_;
};

#endif // PENALTYMERITFUNCTION_H
