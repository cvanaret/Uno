#ifndef UNO_L1MERITFUNCTION_H
#define UNO_L1MERITFUNCTION_H

#include "GlobalizationStrategy.hpp"
#include "tools/Options.hpp"

class l1MeritFunction : public GlobalizationStrategy {
public:
   explicit l1MeritFunction(const Options& options);

   void initialize(Statistics& statistics, const Iterate& first_iterate) override;
   bool check_acceptance(Statistics& statistics, const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress, double
      objective_multiplier, double predicted_reduction) override;
   void reset() override;
   void notify(Iterate& current_iterate) override;

private:
   const double decrease_fraction;
};

#endif // UNO_L1MERITFUNCTION_H
