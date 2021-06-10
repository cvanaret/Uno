#include <iostream>
#include <cmath>
#include "Uno.hpp"
#include "FilterStrategy.hpp"
#include "Vector.hpp"

FilterStrategy::FilterStrategy(ConstraintRelaxationStrategy& feasibility_strategy, Subproblem& subproblem, FilterStrategyParameters&
strategy_parameters, const std::map<std::string, std::string>& options) :
   GlobalizationStrategy(feasibility_strategy, subproblem), filter_optimality(FilterFactory::create(options)),
      filter_restoration(FilterFactory::create(options)), current_phase_(OPTIMALITY), parameters_(strategy_parameters) {
}

Iterate FilterStrategy::initialize(Statistics& statistics, Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
   statistics.add_column("phase", Statistics::int_width, 4);
   /* initialize the subproblem */
   Iterate first_iterate = this->subproblem.evaluate_initial_point(problem, x, multipliers);
   this->subproblem.compute_residuals(problem, first_iterate, first_iterate.multipliers, 1.);
   // preallocate trial_iterate
   this->trial_primals_.resize(first_iterate.x.size());

   /* set the filter upper bound */
   double upper_bound = std::max(this->parameters_.ubd, this->parameters_.fact * first_iterate.progress.feasibility);
   this->filter_optimality->upper_bound = upper_bound;
   this->filter_restoration->upper_bound = upper_bound;
   return first_iterate;
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * precondition: feasible step
 * */
std::optional<Iterate>
FilterStrategy::check_acceptance(Statistics& statistics, Problem& problem, Iterate& current_iterate, Direction& direction,
      double step_length) {
   /* check if subproblem definition changed */
   if (this->subproblem.subproblem_definition_changed) {
      this->filter_optimality->reset();
      this->subproblem.subproblem_definition_changed = false;
      this->subproblem.compute_optimality_measures(problem, current_iterate);
   }

   /* assemble trial point: TODO do not reevaluate if ||d|| = 0 */
   add_vectors(current_iterate.x, direction.x, step_length, this->trial_primals_);
   Iterate trial_iterate(this->trial_primals_, direction.multipliers);
   double step_norm = step_length * direction.norm;
   this->subproblem.compute_optimality_measures(problem, trial_iterate);

   bool accept = false;
   /* check zero step */
   if (step_norm == 0.) {
      accept = true;
   }
   else {
      /* possibly switch phase */
      this->switch_phase_(problem, direction, current_iterate, trial_iterate);

      /* if RESTORATION phase, compute (in)feasibility measures of trial point */
      if (this->current_phase_ == RESTORATION) {
         this->subproblem.compute_infeasibility_measures(problem, trial_iterate, direction);
      }

      DEBUG << "Current: η = " << current_iterate.progress.feasibility << ", ω = " << current_iterate.progress.objective << "\n";
      DEBUG << "Trial: η = " << trial_iterate.progress.feasibility << ", ω = " << trial_iterate.progress.objective << "\n";

      /* check acceptance */
      Filter& filter = (this->current_phase_ == OPTIMALITY) ? *(this->filter_optimality) : *(this->filter_restoration);
      bool acceptable = filter.accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective);
      if (acceptable) {
         // check acceptance wrt current x (h,f)
         acceptable = filter.improves_current_iterate(current_iterate.progress.feasibility, current_iterate.progress.objective,
               trial_iterate.progress.feasibility, trial_iterate.progress.objective);
         if (acceptable) {
            double predicted_reduction = direction.predicted_reduction(step_length);
            double actual_reduction =
                  filter.compute_actual_reduction(current_iterate.progress.objective, current_iterate.progress.feasibility,
                        trial_iterate.progress.objective);
            DEBUG << "Predicted reduction: " << predicted_reduction << ", actual: " << actual_reduction << "\n\n";

            /* switching condition */
            if (predicted_reduction < this->parameters_.Delta * std::pow(current_iterate.progress.feasibility, 2)) {
               filter.add(current_iterate.progress.feasibility, current_iterate.progress.objective);
               accept = true;
            }
               /* Armijo sufficient decrease condition: predicted_reduction should be positive */
            else if (actual_reduction >= this->parameters_.Sigma * step_length * std::max(0., predicted_reduction - 1e-9)) {
               accept = true;
            }
         }
      }
   }

   /* correct multipliers for infeasibility problem */
   if (accept) {
      statistics.add_statistic("phase", (int) direction.phase);
      if (direction.phase == RESTORATION) {
         this->update_restoration_multipliers_(trial_iterate, direction.constraint_partition);
      }
      return trial_iterate;
   }
   else {
      return std::nullopt;
   }
}

void FilterStrategy::switch_phase_(Problem& problem, Direction& direction, Iterate& current_iterate, Iterate& trial_iterate) {
   /* find out if transition of one phase to the other */
   if (this->current_phase_ == OPTIMALITY) {
      if (direction.phase == RESTORATION) {
         /* infeasible QP: go from phase II (optimality) to I (restoration) */
         DEBUG << "Switching from optimality to restoration phase\n";
         this->current_phase_ = RESTORATION;
         /* add [h,f] (c/s violation) to filter, entering restoration */
         this->filter_optimality->add(current_iterate.progress.feasibility, current_iterate.progress.objective);

         /* re-initialize the restoration filter */
         this->filter_restoration->reset();
         this->filter_restoration->upper_bound = this->filter_optimality->upper_bound;
         this->subproblem.compute_infeasibility_measures(problem, current_iterate, direction);
         this->filter_restoration->add(current_iterate.progress.feasibility, current_iterate.progress.objective);
         //current_iterate.is_hessian_computed = false;
      }
   }
      /* check whether we can switch from phase I (restoration) to II (optimality) */
   else if (direction.phase == OPTIMALITY &&
            this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective)) {
      DEBUG << "Switching from restoration to optimality phase\n";
      this->current_phase_ = OPTIMALITY;
      this->subproblem.compute_optimality_measures(problem, current_iterate);
   }
}

void FilterStrategy::update_restoration_multipliers_(Iterate& trial_iterate, ConstraintPartition& constraint_partition) {
   for (int j: constraint_partition.infeasible) {
      if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
         trial_iterate.multipliers.constraints[j] = -1.;
      }
      else {
         trial_iterate.multipliers.constraints[j] = 1.;
      }
   }
}
