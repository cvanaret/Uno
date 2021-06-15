#include <cassert>
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

void FeasibilityRestoration::generate_subproblem(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
   // simply generate the subproblem
   this->subproblem.generate(problem, current_iterate, objective_multiplier, trust_region_radius);
}

Direction FeasibilityRestoration::compute_feasible_direction(const Problem& problem, Iterate& current_iterate) {
   // solve the original subproblem
   Direction direction = this->subproblem.compute_direction(problem, current_iterate);

   if (direction.status != INFEASIBLE) {
      direction.objective_multiplier = problem.objective_sign;
   }
   else {
      ConstraintPartition constraint_partition = direction.constraint_partition;
      // infeasible subproblem: form the feasibility problem
      this->form_feasibility_problem(problem, current_iterate, constraint_partition);
      // solve the feasibility subproblem
      direction = this->subproblem.compute_direction(problem, current_iterate);
      direction.objective_multiplier = 0.;
      direction.constraint_partition = constraint_partition;
      direction.is_relaxed = true;
   }
   return direction;
}

double FeasibilityRestoration::compute_predicted_reduction(const Problem& /*problem*/, Iterate& /*current_iterate*/, Direction& direction,
      double step_length) {
   // the predicted reduction is simply that of the subproblem (the objective multiplier was set accordingly)
   return direction.predicted_reduction(step_length);
}

void FeasibilityRestoration::form_feasibility_problem(const Problem& problem, const Iterate& current_iterate, const ConstraintPartition&
constraint_partition) {
   // set the multipliers of the violated constraints
   this->set_restoration_multipliers(problem, constraint_partition);
   // compute the objective gradient and (possibly) Hessian
   this->subproblem.update_objective_multiplier(problem, current_iterate, 0.);
   this->subproblem.compute_feasibility_linear_objective(current_iterate, constraint_partition);
   this->subproblem.generate_feasibility_bounds(problem, current_iterate.constraints, constraint_partition);
}

Direction FeasibilityRestoration::solve_feasibility_problem(const Problem& problem, Iterate& current_iterate, Direction& direction) {
   assert(this->current_phase == OPTIMALITY && "FeasibilityRestoration is already in the feasibility restoration phase");
   this->form_feasibility_problem(problem, current_iterate, direction.constraint_partition);
   return this->subproblem.compute_direction(problem, current_iterate);
}

bool FeasibilityRestoration::is_acceptable(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Iterate& trial_iterate,
      Direction& direction, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem.subproblem_definition_changed) {
      this->phase_2_strategy->reset();
      this->subproblem.subproblem_definition_changed = false;
      //this->subproblem.compute_optimality_measures(problem, current_iterate);
   }
   double step_norm = step_length * direction.norm;

   this->subproblem.compute_progress_measures(problem, trial_iterate);
   bool accept = false;
   if (step_norm == 0.) {
      accept = true;
   }
   else {
      //ProgressMeasures current_progress_tmp{};
      //ProgressMeasures trial_progress_tmp{};

      /* possibly go from 1 (restoration) to phase 2 (optimality) */
      if (!direction.is_relaxed && this->current_phase == FEASIBILITY_RESTORATION) {
         // TODO && this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective))
         this->current_phase = OPTIMALITY;
         DEBUG << "Switching from restoration to optimality phase\n";
         this->subproblem.compute_progress_measures(problem, current_iterate);
      }
      /* possibly go from phase 2 (optimality) to 1 (restoration) */
      else if (direction.is_relaxed && this->current_phase == OPTIMALITY) {
         this->current_phase = FEASIBILITY_RESTORATION;
         DEBUG << "Switching from optimality to restoration phase\n";
         this->phase_2_strategy->notify(current_iterate);
         this->phase_1_strategy->reset();
         this->compute_infeasibility_measures(problem, current_iterate, direction.constraint_partition);
         this->phase_1_strategy->notify(current_iterate);
      }

      if (this->current_phase == FEASIBILITY_RESTORATION) {
         // if restoration phase, recompute progress measures of trial point
         this->compute_infeasibility_measures(problem, trial_iterate, direction.constraint_partition);
      }

      // evaluate the predicted reduction
      double predicted_reduction = direction.predicted_reduction(step_length);
      // pick the current strategy
      GlobalizationStrategy& strategy = (this->current_phase == OPTIMALITY) ? *this->phase_2_strategy : *this->phase_1_strategy;
      // invoke the globalization strategy for acceptance
      accept = strategy.check_acceptance(statistics, current_iterate.progress, trial_iterate.progress,
            direction.objective_multiplier, predicted_reduction);
   }

   if (accept) {
      statistics.add_statistic("phase", (int) direction.is_relaxed ? FEASIBILITY_RESTORATION : OPTIMALITY);
      if (direction.is_relaxed) {
         /* correct multipliers for infeasibility problem */
         FeasibilityRestoration::update_restoration_multipliers(trial_iterate, direction.constraint_partition);
      }
      // compute the residuals
      trial_iterate.compute_objective(problem);
      this->subproblem.compute_residuals(problem, trial_iterate, trial_iterate.multipliers, direction.objective_multiplier);

   }
   return accept;
}

void FeasibilityRestoration::set_restoration_multipliers(const Problem& problem, const ConstraintPartition& constraint_partition) {
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
         this->subproblem.constraints_multipliers[j] = 1.;
      }
      else if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
         this->subproblem.constraints_multipliers[j] = -1.;
      }
      // otherwise, leave the multiplier as it is
   }
}

void FeasibilityRestoration::compute_infeasibility_measures(const Problem& problem, Iterate& iterate, const ConstraintPartition& constraint_partition) {
   iterate.compute_constraints(problem);
   // feasibility measure: residual of all constraints
   double feasibility = problem.compute_constraint_residual(iterate.constraints, L1_NORM);
   // optimality measure: residual of linearly infeasible constraints
   double objective = problem.compute_constraint_violation(iterate.constraints, constraint_partition.infeasible, L1_NORM);
   iterate.progress = {feasibility, objective};
}

void FeasibilityRestoration::update_restoration_multipliers(Iterate& trial_iterate, const ConstraintPartition& constraint_partition) {
   for (int j: constraint_partition.infeasible) {
      if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
         trial_iterate.multipliers.constraints[j] = -1.;
      }
      else {
         trial_iterate.multipliers.constraints[j] = 1.;
      }
   }
}