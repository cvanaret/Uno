#include <iostream>
#include <cmath>
#include "FilterStrategy.hpp"

FilterStrategy::FilterStrategy(FilterStrategyParameters& strategy_parameters, const Options& options) :
      GlobalizationStrategy(), filter(FilterFactory::create(options)), initial_filter_upper_bound(INFINITY),
      parameters(strategy_parameters) {
}

void FilterStrategy::initialize(Statistics& /*statistics*/, const Iterate& first_iterate) {
   /* set the filter upper bound */
   double upper_bound = std::max(this->parameters.ubd, this->parameters.fact * first_iterate.progress.infeasibility);
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

         /* switching condition violated: predicted reduction is not promising */
         if (!FilterStrategy::switching_condition(predicted_reduction, current_progress.infeasibility, this->parameters.Delta)) {
            filter->add(current_progress.infeasibility, current_progress.objective);
            DEBUG << "Trial iterate was accepted by violating switching condition\n";
            DEBUG << "Current iterate was added to the filter\n";
            accept = true;
         }
         /* Armijo sufficient decrease condition: predicted_reduction should be positive */
         else if (FilterStrategy::armijo_condition(predicted_reduction, actual_reduction, this->parameters.decrease_fraction)) {
            DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
            accept = true;
         }
         else { // switching condition holds, but not Armijo condition
            filter->add(current_progress.infeasibility, current_progress.objective);
            DEBUG << "Armijo condition not satisfied\n";
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

bool FilterStrategy::switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) {
   return predicted_reduction > switching_fraction * std::pow(current_infeasibility, 2);
}

bool FilterStrategy::armijo_condition(double predicted_reduction, double actual_reduction, double decrease_fraction) {
   return actual_reduction >= decrease_fraction * std::max(0., predicted_reduction - 1e-9);
}