// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "FilterStrategy.hpp"
#include "tools/Logger.hpp"

FilterStrategy::FilterStrategy(const Options& options) :
      GlobalizationStrategy(options),
      filter(FilterFactory::create(options)),
      parameters({
         options.get_double("filter_delta"),
         options.get_double("filter_ubd"),
         options.get_double("filter_fact"),
         options.get_double("filter_switching_infeasibility_exponent")
      }) {
}

void FilterStrategy::initialize(const Iterate& first_iterate) {
   // set the filter upper bound
   double upper_bound = std::max(this->parameters.upper_bound, this->parameters.infeasibility_fraction * first_iterate.nonlinear_progress.infeasibility);
   this->filter->upper_bound = upper_bound;
   this->initial_filter_upper_bound = upper_bound;
}

void FilterStrategy::reset() {
   // re-initialize the restoration filter
   this->filter->reset();
   // TODO: we should set the ub of the optimality filter. But now, our 2 filters live independently...
   this->filter->upper_bound = this->initial_filter_upper_bound;
}

void FilterStrategy::notify(Iterate& current_iterate) {
   this->filter->add(current_iterate.nonlinear_progress.infeasibility, current_iterate.nonlinear_progress.optimality);
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * filter methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool FilterStrategy::is_acceptable(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress,
      double /*objective_multiplier*/, const ProgressMeasures& predicted_reduction) {
   GlobalizationStrategy::check_finiteness(current_progress);
   GlobalizationStrategy::check_finiteness(trial_progress);

   // include the predicted optimality reduction, but ignore the predicted infeasibility reduction
   const double unconstrained_predicted_reduction = predicted_reduction.optimality;
   DEBUG << "Current: η = " << current_progress.infeasibility << ", ω = " << current_progress.optimality << '\n';
   DEBUG << "Trial:   η = " << trial_progress.infeasibility << ", ω = " << trial_progress.optimality << '\n';
   DEBUG << "Unconstrained predicted reduction: " << unconstrained_predicted_reduction << '\n';

   bool accept = false;
   // check acceptance
   bool acceptable = filter->accept(trial_progress.infeasibility, trial_progress.optimality);
   if (acceptable) {
      DEBUG << "Filter acceptable\n";
      // check acceptance wrt current x (h,f)
      acceptable = filter->improves_current_iterate(current_progress.infeasibility, current_progress.optimality, trial_progress.infeasibility,
            trial_progress.optimality);
      if (acceptable) {
         const double actual_reduction = filter->compute_actual_reduction(current_progress.optimality, current_progress.infeasibility,
               trial_progress.optimality);
         DEBUG << "Filter acceptable wrt current point\n";
         DEBUG << "Actual reduction: " << actual_reduction << '\n';
         DEBUG << *this->filter << '\n';

         // switching condition violated: predicted reduction is not promising
         if (!this->switching_condition(unconstrained_predicted_reduction, current_progress.infeasibility, this->parameters.delta)) {
            filter->add(current_progress.infeasibility, current_progress.optimality);
            DEBUG << "Trial iterate was accepted by violating switching condition\n";
            DEBUG << "Current iterate was added to the filter\n";
            accept = true;
         }
         // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
         else if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction)) {
            DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
            accept = true;
         }
         else { // switching condition holds, but not Armijo condition
            filter->add(current_progress.infeasibility, current_progress.optimality);
            DEBUG << "Armijo condition not satisfied, trial iterate rejected\n";
            DEBUG << "Current iterate was added to the filter\n";
         }
      }
      else {
         DEBUG << "Not filter acceptable wrt current point\n";
      }
   }
   else {
      DEBUG << "Not filter acceptable\n";
   }
   DEBUG << '\n';
   return accept;
}

bool FilterStrategy::switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const {
   return predicted_reduction > switching_fraction * std::pow(current_infeasibility, this->parameters.switching_infeasibility_exponent);
}
