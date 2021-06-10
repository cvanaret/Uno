#include "FeasibilityRestoration.hpp"
#include "GlobalizationStrategyFactory.hpp"

FeasibilityRestoration::FeasibilityRestoration(const std::string& constraint_relaxation_strategy, Subproblem& subproblem,
   const std::map<std::string, std::string>& options) : ConstraintRelaxationStrategy(subproblem),
   /* create the globalization strategy */
   phase_1_strategy(GlobalizationStrategyFactory::create(constraint_relaxation_strategy, options)),
   phase_2_strategy(GlobalizationStrategyFactory::create(constraint_relaxation_strategy, options)),
   current_phase(OPTIMALITY) {
}

Iterate FeasibilityRestoration::initialize(Statistics& statistics, Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
   statistics.add_column("phase", Statistics::int_width, 4);

   /* initialize the subproblem */
   Iterate first_iterate = this->subproblem.evaluate_initial_point(problem, x, multipliers);
   this->subproblem.compute_residuals(problem, first_iterate, first_iterate.multipliers, 1.);

   // preallocate trial_iterate
   this->trial_primals_.resize(first_iterate.x.size());

   this->phase_1_strategy->initialize(statistics, first_iterate);
   this->phase_2_strategy->initialize(statistics, first_iterate);
   return first_iterate;
}

std::vector<Direction> FeasibilityRestoration::compute_feasible_directions(Problem& problem, Iterate& current_iterate, double
trust_region_radius) {
   std::vector<Direction> directions = this->subproblem.compute_directions(problem, current_iterate, trust_region_radius);
   // TODO handle multiple directions
   Direction& direction = directions[0];

   if (direction.status != INFEASIBLE) {
      direction.phase = OPTIMALITY;
      direction.objective_multiplier = problem.objective_sign;
   }
   else {
      /* infeasible subproblem: switch to restoration phase */
      directions = this->subproblem.restore_feasibility(problem, current_iterate, direction, trust_region_radius);
   }
   return directions;
}

std::optional<Iterate> FeasibilityRestoration::check_acceptance(Statistics& statistics, Problem& problem, Iterate& current_iterate, Direction& direction,
      double step_length) {
   // compute the trial iterate TODO do not reevaluate if ||d|| = 0
   add_vectors(current_iterate.x, direction.x, step_length, this->trial_primals_);
   Iterate trial_iterate(this->trial_primals_, direction.multipliers);
   double step_norm = step_length * direction.norm;
   this->subproblem.compute_optimality_measures(problem, trial_iterate);

   // check if subproblem definition changed
   if (this->subproblem.subproblem_definition_changed) {
      this->phase_2_strategy->reset();
      this->subproblem.subproblem_definition_changed = false;
      this->subproblem.compute_optimality_measures(problem, current_iterate);
   }

   bool accept = false;
   if (step_norm == 0.) {
      accept = true;
   }
   else {
      this->switch_phase_(problem, direction, current_iterate, trial_iterate);

      // evaluate the predicted reduction
      double predicted_reduction = direction.predicted_reduction(step_length);
      if (this->current_phase == OPTIMALITY) {
         accept = this->phase_2_strategy->check_acceptance(statistics, current_iterate.progress, trial_iterate.progress, direction, predicted_reduction);
      }
      else {
         /* if RESTORATION phase, compute (in)feasibility measures of trial point */
         this->subproblem.compute_infeasibility_measures(problem, trial_iterate, direction);
         accept = this->phase_1_strategy->check_acceptance(statistics, current_iterate.progress, trial_iterate.progress, direction, predicted_reduction);
      }
   }

   if (accept) {
      statistics.add_statistic("phase", (int) direction.phase);
      if (direction.phase == FEASIBILITY_RESTORATION) {
         /* correct multipliers for infeasibility problem */
         FeasibilityRestoration::update_restoration_multipliers_(trial_iterate, direction.constraint_partition);
      }
      return trial_iterate;
   }
   else {
      return std::nullopt;
   }
}

void FeasibilityRestoration::switch_phase_(Problem& problem, Direction& direction, Iterate& current_iterate, Iterate& /*trial_iterate*/) {
   /* find out if transition of one phase to the other */
   if (this->current_phase == OPTIMALITY) {
      if (direction.phase == FEASIBILITY_RESTORATION) {
         /* infeasible subproblem: go from phase II (optimality) to I (restoration) */
         DEBUG << "Switching from optimality to restoration phase\n";
         this->current_phase = FEASIBILITY_RESTORATION;

         this->phase_2_strategy->notify(current_iterate);
         this->phase_1_strategy->reset();
         this->subproblem.compute_infeasibility_measures(problem, current_iterate, direction);
         this->phase_1_strategy->notify(current_iterate);
      }
   }
      /* check whether we can switch from phase I (restoration) to II (optimality) */
   else if (direction.phase == OPTIMALITY) { // TODO && this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress
   // .objective)) {
      DEBUG << "Switching from restoration to optimality phase\n";
      this->current_phase = OPTIMALITY;
      this->subproblem.compute_optimality_measures(problem, current_iterate);
   }
}

double FeasibilityRestoration::compute_predicted_reduction(Problem& /*problem*/, Iterate& /*current_iterate*/, Direction& direction,
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