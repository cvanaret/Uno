#include "WaechterFilterStrategy.hpp"

WaechterFilterStrategy::WaechterFilterStrategy(const Options& options): FilterStrategy(options) {
}

void WaechterFilterStrategy::initialize(const Iterate& first_iterate) {
   this->initial_infeasibility = first_iterate.residuals.infeasibility;
   FilterStrategy::initialize(first_iterate);
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * filter methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool WaechterFilterStrategy::is_iterate_acceptable(const ProgressMeasures& current_progress_measures, const ProgressMeasures& trial_progress_measures,
      const PredictedReduction& predicted_reduction, double /*objective_multiplier*/) {
   const double current_optimality_measure = current_progress_measures.scaled_optimality(1.) + current_progress_measures.unscaled_optimality;
   const double trial_optimality_measure = trial_progress_measures.scaled_optimality(1.) + trial_progress_measures.unscaled_optimality;
   // unconstrained predicted reduction:
   // - ignore the predicted infeasibility reduction
   // - scale the scaled optimality measure with 1
   const double unconstrained_predicted_reduction = predicted_reduction.scaled_optimality(1.) + predicted_reduction.unscaled_optimality;
   DEBUG << "Current: η = " << current_progress_measures.infeasibility << ", ω = " << current_optimality_measure << '\n';
   DEBUG << "Trial:   η = " << trial_progress_measures.infeasibility << ", ω = " << trial_optimality_measure << '\n';
   DEBUG << "Unconstrained predicted reduction: " << predicted_reduction.scaled_optimality(1.) << " + " << predicted_reduction.unscaled_optimality <<
         " = " <<  unconstrained_predicted_reduction << '\n';

   GlobalizationStrategy::check_finiteness(current_progress_measures, 1.);
   GlobalizationStrategy::check_finiteness(trial_progress_measures, 1.);
   DEBUG << *this->filter << '\n';

   bool accept = false;
   // check acceptance
   const bool filter_acceptable = this->filter->accept(trial_progress_measures.infeasibility, trial_optimality_measure);
   if (filter_acceptable) {
      DEBUG << "Filter acceptable\n";
      // compute actual reduction (and protect against roundoff errors)
      static double machine_epsilon = std::numeric_limits<double>::epsilon();
      double actual_reduction = this->filter->compute_actual_reduction(current_optimality_measure, current_progress_measures.infeasibility,
            trial_optimality_measure) + 10. * machine_epsilon * std::abs(current_optimality_measure);
      DEBUG << "Actual reduction: " << actual_reduction << '\n';

      const bool small_infeasibility = current_progress_measures.infeasibility <= 1e-4*std::max(1., this->initial_infeasibility);
      const bool switching = (0. < unconstrained_predicted_reduction) && this->switching_condition(unconstrained_predicted_reduction,
            current_progress_measures.infeasibility, this->parameters.delta);
      if (small_infeasibility && switching) {
         DEBUG << "Switching condition satisfied\n";
         // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
         if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction)) {
            DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
            accept = true;
         }
         else {
            DEBUG << "Armijo condition not satisfied\n";
            this->filter->add(current_progress_measures.infeasibility, current_optimality_measure);
         }
      }
      if (not small_infeasibility or not switching) {
         DEBUG << "Switching condition violated\n";
         if (this->filter->improves_current_iterate(current_progress_measures.infeasibility, current_optimality_measure,
               trial_progress_measures.infeasibility, trial_optimality_measure)) {
            DEBUG << "Acceptable wrt current point\n";
            accept = true;
         }
         else {
            DEBUG << "Not acceptable wrt current point\n";
         }
         this->filter->add(current_progress_measures.infeasibility, current_optimality_measure);
      }
   }
   else {
      DEBUG << "Not filter acceptable\n";
   }
   DEBUG << '\n';
   return accept;
}