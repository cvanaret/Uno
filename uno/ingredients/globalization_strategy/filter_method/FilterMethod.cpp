// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FilterMethod.hpp"
#include "filter/FilterFactory.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Options.hpp"
#include "tools/Statistics.hpp"

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

/* check acceptability of step(s) (filter & sufficient reduction)
 * filter methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool FilterMethod::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
      const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) {
   const bool solving_feasibility_problem = (objective_multiplier == 0.);
   if (solving_feasibility_problem) {
      return this->is_feasibility_iterate_acceptable(statistics, current_progress, trial_progress, predicted_reduction);
   }
   else {
      return this->is_regular_iterate_acceptable(statistics, current_progress, trial_progress, predicted_reduction);
   }
}

void FilterMethod::reset() {
   this->filter->reset();
}

void FilterMethod::register_current_progress(const ProgressMeasures& current_progress) {
   const double current_objective_measure = FilterMethod::unconstrained_merit_function(current_progress);
   this->filter->add(current_progress.infeasibility, current_objective_measure);
}

// solving the feasibility problem = working on infeasibility only (no filter acceptability test)
bool FilterMethod::is_feasibility_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
      const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction) const {
   DEBUG << "Current infeasibility = " << current_progress.infeasibility << '\n';
   DEBUG << "Trial   infeasibility = " << trial_progress.infeasibility << '\n';
   DEBUG << "Predicted reduction = " << predicted_reduction.infeasibility << '\n';
   bool accept = false;
   if (this->armijo_sufficient_decrease(predicted_reduction.infeasibility, current_progress.infeasibility - trial_progress.infeasibility)) {
      DEBUG << "Trial iterate (h-type) was accepted by satisfying the Armijo condition\n";
      accept = true;
   }
   else {
      DEBUG << "Trial iterate (h-type) was rejected by violating the Armijo condition\n";
   }
   Iterate::number_eval_objective--;
   statistics.set("status", std::string(accept ? "accepted" : "rejected") + " (h-type Armijo)");
   return accept;
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
