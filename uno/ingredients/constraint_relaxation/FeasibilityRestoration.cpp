#include <cassert>
#include "FeasibilityRestoration.hpp"
#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

FeasibilityRestoration::FeasibilityRestoration(const Problem& problem, Subproblem& subproblem, const Options& options) :
      ConstraintRelaxationStrategy(problem, subproblem),
      // create the globalization strategy
      phase_1_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)),
      phase_2_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)) {
}

void FeasibilityRestoration::initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) {
   statistics.add_column("phase", Statistics::int_width, 4);

   // initialize the subproblem
   this->subproblem.initialize(statistics, problem, first_iterate);
   this->subproblem.compute_errors(problem, first_iterate, problem.objective_sign);

   this->phase_1_strategy->initialize(statistics, first_iterate);
   this->phase_2_strategy->initialize(statistics, first_iterate);
}

void FeasibilityRestoration::create_current_subproblem(const Problem& problem, Iterate& current_iterate, double trust_region_radius) {
   // simply generate the subproblem
   this->subproblem.create_current_subproblem(problem, current_iterate, problem.objective_sign, trust_region_radius);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   // solve the original subproblem
   Direction direction = this->subproblem.solve(statistics, problem, current_iterate);
   direction.objective_multiplier = problem.objective_sign;
   DEBUG << "\n" << direction;

   // infeasible subproblem: form the feasibility problem
   if (direction.status == INFEASIBLE) {
      // try to minimize the constraint violation by solving the feasibility subproblem
      direction = this->solve_feasibility_problem(statistics, problem, current_iterate, direction);
      DEBUG << "\n" << direction;
   }
   return direction;
}

double FeasibilityRestoration::compute_predicted_reduction(const Problem& /*problem*/, Iterate& /*current_iterate*/, const Direction& /*direction*/,
      PredictedReductionModel& predicted_reduction_model, double step_length) {
   // the predicted reduction is simply that of the subproblem (the objective multiplier was set accordingly)
   return predicted_reduction_model.evaluate(step_length);
}

void FeasibilityRestoration::form_feasibility_problem(const Problem& problem, Iterate& current_iterate, const std::vector<double>&
phase_2_primal_direction, const ConstraintPartition& constraint_partition) {
   assert(!constraint_partition.infeasible.empty() && "The direction is infeasible but no constraint is infeasible");

   // set the multipliers of the violated constraints
   FeasibilityRestoration::set_restoration_multipliers(this->subproblem.constraints_multipliers, constraint_partition);

   // build the local model of the objective
   this->subproblem.build_objective_model(problem, current_iterate, 0.);
   this->subproblem.compute_feasibility_linear_objective(current_iterate, constraint_partition);

   // update the bounds of the constraints
   this->subproblem.generate_feasibility_bounds(problem, current_iterate.constraints, constraint_partition);

   // start from the phase-2 solution
   this->subproblem.set_initial_point(phase_2_primal_direction);
}

void FeasibilityRestoration::form_feasibility_problem(const Problem& problem, Iterate& current_iterate, const std::vector<double>&
phase_2_primal_direction) {
   assert(false && "form_feasibility_problem without constraint partition: must be implemented");

   // compute the objective model with a zero objective multiplier
   this->subproblem.build_objective_model(problem, current_iterate, 0.);

   // add elastic variables to relax the problem
   this->add_elastic_variables_to_subproblem();

   // start from the phase-2 solution
   this->subproblem.set_initial_point(phase_2_primal_direction);
}

Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, const Problem& problem, Iterate& current_iterate,
      const Direction& phase_2_direction) {
   // form the feasibility problem (with or without constraint partition)
   if (phase_2_direction.constraint_partition.has_value()) {
      ConstraintPartition constraint_partition = phase_2_direction.constraint_partition.value();
      this->form_feasibility_problem(problem, current_iterate, phase_2_direction.x, constraint_partition);
   }
   else {
      this->form_feasibility_problem(problem, current_iterate, phase_2_direction.x);
   }

   // solve the feasibility subproblem
   DEBUG << "\nSolving the feasibility subproblem\n";
   Direction feasibility_direction = this->subproblem.solve(statistics, problem, current_iterate);
   feasibility_direction.objective_multiplier = 0.;
   feasibility_direction.is_relaxed = true;
   if (phase_2_direction.constraint_partition.has_value()) {
      ConstraintPartition constraint_partition = phase_2_direction.constraint_partition.value();
      feasibility_direction.constraint_partition = constraint_partition;
   }
   else {
      // TODO remove elastic variables
      std::cout << feasibility_direction << "\n";
      assert(false && "TODO");
   }
   return feasibility_direction;
}

bool FeasibilityRestoration::is_acceptable(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Iterate& trial_iterate,
      const Direction& direction, PredictedReductionModel& predicted_reduction_model, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem.subproblem_definition_changed) {
      this->phase_2_strategy->reset();
      this->subproblem.subproblem_definition_changed = false;
      this->subproblem.compute_progress_measures(problem, current_iterate);
   }

   bool accept = false;
   if (direction.norm == 0.) {
      accept = true;
   }
   else {
      // possibly go from 1 (restoration) to phase 2 (optimality)
      if (!direction.is_relaxed && this->current_phase == FEASIBILITY_RESTORATION) {
         // TODO && this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective))
         this->current_phase = OPTIMALITY;
         DEBUG << "Switching from restoration to optimality phase\n";
         this->subproblem.compute_progress_measures(problem, current_iterate);
      }
         // possibly go from phase 2 (optimality) to 1 (restoration)
      else if (direction.is_relaxed && this->current_phase == OPTIMALITY) {
         this->current_phase = FEASIBILITY_RESTORATION;
         DEBUG << "Switching from optimality to restoration phase\n";
         this->phase_2_strategy->notify(current_iterate);
         this->phase_1_strategy->reset();
         this->compute_infeasibility_measures(problem, current_iterate, direction.constraint_partition);
         this->phase_1_strategy->notify(current_iterate);
      }

      trial_iterate.evaluate_constraints(problem);
      if (this->current_phase == FEASIBILITY_RESTORATION) {
         // if restoration phase, recompute progress measures of trial point
         this->compute_infeasibility_measures(problem, trial_iterate, direction.constraint_partition);
      }
      else {
         this->subproblem.compute_progress_measures(problem, trial_iterate);
      }

      // evaluate the predicted reduction
      const double predicted_reduction = this->compute_predicted_reduction(problem, current_iterate, direction, predicted_reduction_model, step_length);

      // pick the current strategy
      GlobalizationStrategy& strategy = (this->current_phase == OPTIMALITY) ? *this->phase_2_strategy : *this->phase_1_strategy;
      // invoke the globalization strategy for acceptance
      accept = strategy.check_acceptance(statistics, current_iterate.progress, trial_iterate.progress,
            direction.objective_multiplier, predicted_reduction);
   }

   if (accept) {
      statistics.add_statistic("phase", static_cast<int>(direction.is_relaxed ? FEASIBILITY_RESTORATION : OPTIMALITY));
      if (direction.is_relaxed && direction.constraint_partition.has_value()) {
         // correct multipliers for infeasibility problem
         FeasibilityRestoration::set_restoration_multipliers(trial_iterate.multipliers.constraints, direction.constraint_partition.value());
      }
      // compute the errors
      trial_iterate.evaluate_objective(problem);
      this->subproblem.compute_errors(problem, trial_iterate, direction.objective_multiplier);
   }
   return accept;
}

void FeasibilityRestoration::set_restoration_multipliers(std::vector<double>& constraints_multipliers, const ConstraintPartition&
constraint_partition) {
   for (size_t j: constraint_partition.infeasible) {
      if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
         constraints_multipliers[j] = 1.;
      }
      else { // constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER
         constraints_multipliers[j] = -1.;
      }
   }
   // otherwise, leave the multiplier as it is
}

void FeasibilityRestoration::compute_infeasibility_measures(const Problem& problem, Iterate& iterate,
      const std::optional<ConstraintPartition>& constraint_partition) {
   iterate.evaluate_constraints(problem);
   // optimality measure: residual of linearly infeasible constraints
   if (constraint_partition.has_value()) {
      // feasibility measure: residual of all constraints
      double feasibility = problem.compute_constraint_violation(iterate.constraints, L1_NORM);
      double objective = problem.compute_constraint_violation(iterate.constraints, constraint_partition.value().infeasible, L1_NORM);
      iterate.progress = {feasibility, objective};
   }
   else {
      // if no constraint partition is available, simply compute the standard progress measures
      this->subproblem.compute_progress_measures(problem, iterate);
   }
}

size_t FeasibilityRestoration::get_max_number_variables(const Problem& problem) {
   return problem.number_variables;
}