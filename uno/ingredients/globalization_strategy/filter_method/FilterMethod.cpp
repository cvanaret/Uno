// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FilterMethod.hpp"
#include "filter/FilterFactory.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Options.hpp"

FilterMethod::FilterMethod(const Options& options) :
      GlobalizationStrategy(options),
      filter(FilterFactory::create(options)),
      parameters({
         options.get_double("filter_delta"),
         options.get_double("filter_ubd"),
         options.get_double("filter_fact"),
         options.get_double("filter_switching_infeasibility_exponent"),
      }) {
}

void FilterMethod::initialize(Statistics& /*statistics*/, const Iterate& initial_iterate, const Options& /*options*/) {
   // set the filter upper bound
   double upper_bound = std::max(this->parameters.upper_bound, this->parameters.infeasibility_fraction * initial_iterate.progress.infeasibility);
   this->filter->set_infeasibility_upper_bound(upper_bound);
}

void FilterMethod::reset() {
   this->filter->reset();
}

void FilterMethod::register_current_progress(const ProgressMeasures& current_progress) {
   const double current_objective_measure = FilterMethod::unconstrained_merit_function(current_progress);
   this->filter->add(current_progress.infeasibility, current_objective_measure);
}

double FilterMethod::unconstrained_merit_function(const ProgressMeasures& progress) {
   return progress.objective(1.) + progress.auxiliary;
}

double FilterMethod::compute_actual_objective_reduction(double current_objective_measure, double current_infeasibility, double trial_objective_measure) {
   double actual_reduction = this->filter->compute_actual_objective_reduction(current_objective_measure, current_infeasibility, trial_objective_measure);
   if (this->protect_actual_reduction_against_roundoff) {
      static double machine_epsilon = std::numeric_limits<double>::epsilon();
      actual_reduction += 10. * machine_epsilon * std::abs(current_objective_measure);
   }
   return actual_reduction;
}

bool FilterMethod::switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const {
   return predicted_reduction > switching_fraction * std::pow(current_infeasibility, this->parameters.switching_infeasibility_exponent);
}
