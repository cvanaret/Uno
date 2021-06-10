#include <iostream>
#include <cmath>
#include "Uno.hpp"
#include "FilterStrategy.hpp"
#include "Vector.hpp"

FilterStrategy::FilterStrategy(Subproblem& subproblem, FilterStrategyParameters& strategy_parameters, const std::map<std::string, std::string>&
      options) :
   GlobalizationStrategy(subproblem), filter(FilterFactory::create(options)), initial_filter_upper_bound(INFINITY),
   parameters_(strategy_parameters) {
}

void FilterStrategy::initialize(Statistics& /*statistics*/, const Iterate& first_iterate) {
   /* set the filter upper bound */
   double upper_bound = std::max(this->parameters_.ubd, this->parameters_.fact * first_iterate.progress.feasibility);
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
   this->filter->add(current_iterate.progress.feasibility, current_iterate.progress.objective);
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * precondition: feasible step
 * */
bool FilterStrategy::check_acceptance(Statistics& /*statistics*/, Iterate& current_iterate, Iterate& trial_iterate, Direction& direction,
      double step_length) {
   bool accept = false;

   DEBUG << "Current: η = " << current_iterate.progress.feasibility << ", ω = " << current_iterate.progress.objective << "\n";
   DEBUG << "Trial: η = " << trial_iterate.progress.feasibility << ", ω = " << trial_iterate.progress.objective << "\n";

   /* check acceptance */
   bool acceptable = filter->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective);
   if (acceptable) {
      // check acceptance wrt current x (h,f)
      acceptable = filter->improves_current_iterate(current_iterate.progress.feasibility, current_iterate.progress.objective,
            trial_iterate.progress.feasibility, trial_iterate.progress.objective);
      if (acceptable) {
         double predicted_reduction = direction.predicted_reduction(step_length);
         double actual_reduction =
               filter->compute_actual_reduction(current_iterate.progress.objective, current_iterate.progress.feasibility,
                     trial_iterate.progress.objective);
         DEBUG << "Predicted reduction: " << predicted_reduction << ", actual: " << actual_reduction << "\n\n";

         /* switching condition */
         if (predicted_reduction < this->parameters_.Delta * std::pow(current_iterate.progress.feasibility, 2)) {
            filter->add(current_iterate.progress.feasibility, current_iterate.progress.objective);
            accept = true;
         }
            /* Armijo sufficient decrease condition: predicted_reduction should be positive */
         else if (actual_reduction >= this->parameters_.Sigma * step_length * std::max(0., predicted_reduction - 1e-9)) {
            accept = true;
         }
      }
   }
   return accept;
}