#include "FeasibilityRestoration.hpp"
#include "GlobalizationStrategyFactory.hpp"

FeasibilityRestoration::FeasibilityRestoration(const std::string& constraint_relaxation_strategy, Subproblem& subproblem,
   const std::map<std::string, std::string>& options) : ConstraintRelaxationStrategy(subproblem),
   /* create the globalization strategy */
   phase_1_strategy(GlobalizationStrategyFactory::create(constraint_relaxation_strategy, options)),
   phase_2_strategy(GlobalizationStrategyFactory::create(constraint_relaxation_strategy, options)),
   current_phase(OPTIMALITY) {
}

Iterate FeasibilityRestoration::initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
   statistics.add_column("phase", Statistics::int_width, 4);

   /* initialize the subproblem */
   Iterate first_iterate = this->subproblem.evaluate_initial_point(problem, x, multipliers);
   this->subproblem.compute_residuals(problem, first_iterate, first_iterate.multipliers, 1.);

   this->phase_1_strategy->initialize(statistics, first_iterate);
   this->phase_2_strategy->initialize(statistics, first_iterate);
   return first_iterate;
}

Direction FeasibilityRestoration::compute_feasible_direction(const Problem& problem, Iterate& current_iterate, double trust_region_radius) {
   Direction direction = this->subproblem.compute_direction(problem, current_iterate, trust_region_radius);

   if (direction.status != INFEASIBLE) {
      direction.objective_multiplier = problem.objective_sign;

      /* possibly go from 1 (restoration) to phase 2 (optimality) */
      if (this->current_phase == FEASIBILITY_RESTORATION) {
         // TODO && this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective))
         DEBUG << "Switching from restoration to optimality phase\n";
         this->current_phase = OPTIMALITY;
         this->subproblem.compute_optimality_measures(problem, current_iterate);
      }
   }
   else {
      /* infeasible subproblem: enter restoration phase */
      current_iterate.multipliers.objective = 0.;
      direction = this->subproblem.restore_feasibility(problem, current_iterate, direction, trust_region_radius);
      direction.objective_multiplier = 0.;
      direction.is_relaxed = true;

      /* possibly go from phase 2 (optimality) to 1 (restoration) */
      if (this->current_phase == OPTIMALITY) {
         DEBUG << "Switching from optimality to restoration phase\n";
         this->current_phase = FEASIBILITY_RESTORATION;

         this->phase_2_strategy->notify(current_iterate);
         this->phase_1_strategy->reset();
         this->subproblem.compute_infeasibility_measures(problem, current_iterate, direction);
         this->phase_1_strategy->notify(current_iterate);
      }
   }
   return direction;
}

bool FeasibilityRestoration::is_acceptable(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Iterate& trial_iterate,
      Direction& direction, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem.subproblem_definition_changed) {
      this->phase_2_strategy->reset();
      this->subproblem.subproblem_definition_changed = false;
      this->subproblem.compute_optimality_measures(problem, current_iterate);
   }
   double step_norm = step_length * direction.norm;

   bool accept = false;
   if (step_norm == 0.) {
      accept = true;
   }
   else {
      if (this->current_phase != OPTIMALITY) {
         // if restoration phase, recompute progress measures of trial point
         this->subproblem.compute_infeasibility_measures(problem, trial_iterate, direction);
      }

      // evaluate the predicted reduction
      double predicted_reduction = direction.predicted_reduction(step_length);
      GlobalizationStrategy& globalization_strategy = (this->current_phase == OPTIMALITY) ?  *this->phase_2_strategy : *this->phase_1_strategy;
      // invoke the globalization strategy for acceptance
      accept = globalization_strategy.check_acceptance(statistics, current_iterate.progress, trial_iterate.progress,
            direction.objective_multiplier, predicted_reduction);
   }

   if (accept) {
      statistics.add_statistic("phase", (int) direction.is_relaxed ? FEASIBILITY_RESTORATION : OPTIMALITY);
      if (direction.is_relaxed) {
         /* correct multipliers for infeasibility problem */
         FeasibilityRestoration::update_restoration_multipliers_(trial_iterate, direction.constraint_partition);
      }
   }
   return accept;
}

double FeasibilityRestoration::compute_predicted_reduction(const Problem& /*problem*/, Iterate& /*current_iterate*/, Direction& direction,
      double step_length) {
   return direction.predicted_reduction(step_length);
}

void FeasibilityRestoration::update_restoration_multipliers_(Iterate& trial_iterate, const ConstraintPartition& constraint_partition) {
   for (int j: constraint_partition.infeasible) {
      if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
         trial_iterate.multipliers.constraints[j] = -1.;
      }
      else {
         trial_iterate.multipliers.constraints[j] = 1.;
      }
   }
}