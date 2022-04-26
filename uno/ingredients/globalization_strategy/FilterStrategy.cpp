#include <cmath>
#include "FilterStrategy.hpp"
#include "tools/Logger.hpp"

FilterStrategy::FilterStrategy(const Options& options) :
      GlobalizationStrategy(options),
      filter(FilterFactory::create(options)),
      parameters({stod(options.at("filter_delta")),
                  stod(options.at("filter_ubd")),
                  stod(options.at("filter_fact")),
                  stod(options.at("filter_switching_infeasibility_exponent"))}) {
}

void FilterStrategy::initialize(Statistics& /*statistics*/, const Iterate& first_iterate) {
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
   this->filter->add(current_iterate.nonlinear_progress.infeasibility, current_iterate.nonlinear_progress.reformulation_objective);
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * precondition: feasible step
 * */
bool FilterStrategy::is_acceptable(Statistics& /*statistics*/, const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress,
      double /*objective_multiplier*/, double predicted_reduction) {
   GlobalizationStrategy::check_finiteness(current_progress);
   GlobalizationStrategy::check_finiteness(trial_progress);

   DEBUG << "Current: η = " << current_progress.infeasibility << ", ω = " << current_progress.reformulation_objective << '\n';
   DEBUG << "Trial:   η = " << trial_progress.infeasibility << ", ω = " << trial_progress.reformulation_objective << '\n';
   DEBUG << "Predicted reduction: " << predicted_reduction << '\n';

   bool accept = false;
   // check acceptance
   bool acceptable = filter->accept(trial_progress.infeasibility, trial_progress.reformulation_objective);
   if (acceptable) {
      DEBUG << "Filter acceptable\n";
      // check acceptance wrt current x (h,f)
      acceptable = filter->improves_current_iterate(current_progress.infeasibility, current_progress.reformulation_objective, trial_progress.infeasibility, trial_progress.reformulation_objective);
      if (acceptable) {
         DEBUG << "Filter acceptable wrt current point\n";
         const double actual_reduction = filter->compute_actual_reduction(current_progress.reformulation_objective, current_progress.infeasibility, trial_progress
         .reformulation_objective);
         DEBUG << "Actual reduction: " << actual_reduction << '\n';
         DEBUG << *this->filter << '\n';

         // switching condition violated: predicted reduction is not promising
         if (!this->switching_condition(predicted_reduction, current_progress.infeasibility, this->parameters.delta)) {
            filter->add(current_progress.infeasibility, current_progress.reformulation_objective);
            DEBUG << "Trial iterate was accepted by violating switching condition\n";
            DEBUG << "Current iterate was added to the filter\n";
            accept = true;
         }
         // Armijo sufficient decrease condition (predicted_reduction should be positive)
         else if (this->armijo_sufficient_decrease(predicted_reduction, actual_reduction)) {
            DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
            accept = true;
         }
         else { // switching condition holds, but not Armijo condition
            filter->add(current_progress.infeasibility, current_progress.reformulation_objective);
            DEBUG << "Armijo condition not satisfied\n";
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