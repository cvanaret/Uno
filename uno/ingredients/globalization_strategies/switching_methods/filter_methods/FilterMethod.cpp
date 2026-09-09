// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FilterMethod.hpp"
#include "filters/Filter.hpp"
#include "filters/FilterFactory.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"

namespace uno {
   FilterMethod::FilterMethod(const Options& options) :
         SwitchingMethod(options),
         filter(FilterFactory::create(options)),
         parameters({
            options.get_double("filter_ubd"),
            options.get_double("filter_fact"),
         }) {
   }

   FilterMethod::~FilterMethod() = default;

   void FilterMethod::initialize(Statistics& /*statistics*/, const Iterate& initial_iterate) {
      // set the filter upper bound
      const double upper_bound = std::max(this->parameters.upper_bound, this->parameters.infeasibility_factor * initial_iterate.progress.infeasibility);
      this->filter->set_infeasibility_upper_bound(upper_bound);
      this->reset();
   }

   void FilterMethod::reset() {
      this->filter->reset();
   }

   void FilterMethod::avoid_cycling_back_to(const ProgressMeasures& current_progress) {
      const double current_objective_measure = unconstrained_merit_function(current_progress);
      this->filter->add(current_progress.infeasibility, current_objective_measure);
   }

   double FilterMethod::compute_actual_objective_reduction(double current_objective_measure, double current_infeasibility, double trial_objective_measure) const {
      double actual_reduction = this->filter->compute_actual_objective_reduction(current_objective_measure, current_infeasibility, trial_objective_measure);
      if (this->protect_actual_reduction_against_roundoff) {
         static double machine_epsilon = std::numeric_limits<double>::epsilon();
         actual_reduction += this->protected_actual_reduction_macheps_coefficient * machine_epsilon * std::abs(current_objective_measure);
      }
      return actual_reduction;
   }
} // namespace