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

void FilterStrategy::register_current_progress(const ProgressMeasures& current_progress_measures) {
   const double current_optimality_measure = current_progress_measures.scaled_optimality(1.) + current_progress_measures.unscaled_optimality;
   this->filter->add(current_progress_measures.infeasibility, current_optimality_measure);
}

bool FilterStrategy::is_feasibility_iterate_acceptable(double trial_infeasibility_measure) const {
   if (this->filter->is_empty()) {
      return true;
   }
   else { // filter not empty
      // accept if the infeasibility measure improves upon the smallest filter infeasibility
      return (trial_infeasibility_measure < this->filter->get_smallest_infeasibility());
   }
}

bool FilterStrategy::switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const {
   return predicted_reduction > switching_fraction * std::pow(current_infeasibility, this->parameters.switching_infeasibility_exponent);
}
