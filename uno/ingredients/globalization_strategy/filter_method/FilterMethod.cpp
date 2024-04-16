// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FilterMethod.hpp"
#include "filter/FilterFactory.hpp"

FilterMethod::FilterMethod(const Options& options) :
      GlobalizationStrategy(options),
      filter(FilterFactory::create(options)),
      parameters({
         options.get_double("filter_delta"),
         options.get_double("filter_ubd"),
         options.get_double("filter_fact"),
         options.get_double("filter_switching_infeasibility_exponent")
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

void FilterMethod::register_current_progress(const ProgressMeasures& current_progress_measures) {
   const double current_objective_measure = current_progress_measures.objective(1.) + current_progress_measures.auxiliary;
   this->filter->add(current_progress_measures.infeasibility, current_objective_measure);
}

double FilterMethod::get_infeasibility_upper_bound() const {
   return this->filter->get_infeasibility_upper_bound();
}

void FilterMethod::set_infeasibility_upper_bound(double new_upper_bound) {
   this->filter->set_infeasibility_upper_bound(new_upper_bound);
}

bool FilterMethod::switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const {
   return predicted_reduction > switching_fraction * std::pow(current_infeasibility, this->parameters.switching_infeasibility_exponent);
}
