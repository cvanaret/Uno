#include <iostream>
#include <cmath>
#include "Uno.hpp"
#include "FilterStrategy.hpp"
#include "Vector.hpp"

FilterStrategy::FilterStrategy(FilterStrategyParameters& strategy_parameters, const Options& options) :
      GlobalizationStrategy(), filter(FilterFactory::create(options)), initial_filter_upper_bound(INFINITY),
      parameters_(strategy_parameters) {
}

void FilterStrategy::initialize(Statistics& /*statistics*/, const Iterate& first_iterate) {
   /* set the filter upper bound */
   double upper_bound = std::max(this->parameters_.ubd, this->parameters_.fact * first_iterate.progress.infeasibility);
   this->filter->upper_bound = upper_bound;
   this->initial_filter_upper_bound = upper_bound;
}

void FilterStrategy::reset() {
   /* re-initialize the restoration filter */
   this->filter->reset();
   // TODO: we should set the ub of the optimality filter. But now, our 2 filters live independently...
   this->filter->upper_bound = this->initial_filter_upper_bound;
}

void FilterStrategy::notify(Iterate& current_iterate) {
   this->filter->add(current_iterate.progress.infeasibility, current_iterate.progress.objective);
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * precondition: feasible step
 * */
bool FilterStrategy::check_acceptance(Statistics& /*statistics*/, ProgressMeasures& current_progress, ProgressMeasures& trial_progress,
      double /*objective_multiplier*/, double predicted_reduction) {
   DEBUG << "Current: η = " << current_progress.infeasibility << ", ω = " << current_progress.objective << "\n";
   DEBUG << "Trial:   η = " << trial_progress.infeasibility << ", ω = " << trial_progress.objective << "\n";
   DEBUG << "Predicted reduction (should be positive): " << predicted_reduction << "\n";
   DEBUG << *this->filter << "\n";

   bool accept = false;
   /* check acceptance */
   bool acceptable = filter->accept(trial_progress.infeasibility, trial_progress.objective);
   if (acceptable) {
      // check acceptance wrt current x (h,f)
      acceptable = filter->improves_current_iterate(current_progress.infeasibility, current_progress.objective, trial_progress.infeasibility, trial_progress.objective);
      if (acceptable) {
         const double actual_reduction = filter->compute_actual_reduction(current_progress.objective, current_progress.infeasibility, trial_progress
         .objective);
         DEBUG << "Actual reduction: " << actual_reduction << "\n";

         /* reverse switching condition: predicted reduction is not promising */
         if (predicted_reduction < this->parameters_.Delta * std::pow(current_progress.infeasibility, 2)) {
            filter->add(current_progress.infeasibility, current_progress.objective);
            DEBUG << "Trial iterate was accepted by switching condition\n";
            DEBUG << "Current iterate was added to the filter\n";
            accept = true;
         }
         /* Armijo sufficient decrease condition: predicted_reduction should be positive */
         else if (actual_reduction >= this->parameters_.Sigma * std::max(0., predicted_reduction - 1e-9)) {
            DEBUG << "Trial iterate was accepted by Armijo condition\n";
            accept = true;
         }
         else {
            DEBUG << "Armijo condition not satisfied: " << actual_reduction << " < " << this->parameters_.Sigma * std::max(0., predicted_reduction
            - 1e-9) << "\n";
         }
      }
      else {
         DEBUG << "Not filter acceptable wrt current point\n";
      }
   }
   else {
      DEBUG << "Not filter acceptable\n";
   }
   return accept;
}