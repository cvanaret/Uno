#include <cassert>
#include "FeasibilityRestoration.hpp"
#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

FeasibilityRestoration::FeasibilityRestoration(const Problem& problem, const Options& options) :
      ConstraintRelaxationStrategy(problem, options),
      // create the globalization strategies (one for each phase)
      phase_1_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)),
      phase_2_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)) {
}

void FeasibilityRestoration::initialize(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& first_iterate) {
   statistics.add_column("phase", Statistics::int_width, 4);

   // initialize the subproblem
   this->subproblem->initialize(statistics, problem, scaling, first_iterate);

   // compute the progress measures and the residuals of the initial point
   this->subproblem->compute_progress_measures(problem, scaling, first_iterate);
   this->subproblem->compute_residuals(problem, scaling, first_iterate, problem.objective_sign);

   // initialize the globalization strategies
   this->phase_1_strategy->initialize(statistics, first_iterate);
   this->phase_2_strategy->initialize(statistics, first_iterate);
}

void FeasibilityRestoration::create_current_subproblem(const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
      double trust_region_radius) {
   this->subproblem->create_current_subproblem(problem, scaling, current_iterate, problem.objective_sign, trust_region_radius);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, const Problem& problem, const Scaling& scaling,
      Iterate& current_iterate) {
   // solve the original subproblem
   Direction direction = this->subproblem->solve(statistics, problem, current_iterate);
   direction.objective_multiplier = problem.objective_sign;
   DEBUG << direction << "\n";
   assert((direction.status == OPTIMAL || direction.status == INFEASIBLE) && "The subproblem was not solved properly");

   // infeasible subproblem: form the feasibility problem
   if (direction.status == INFEASIBLE) {
      // try to minimize the constraint violation by solving the feasibility subproblem
      direction = this->solve_feasibility_problem(statistics, problem, scaling, current_iterate, direction.x, direction.constraint_partition);
      DEBUG << direction << "\n";
   }
   return direction;
}

Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, const Problem& problem, const Scaling& scaling,
      Iterate& current_iterate, const std::optional<std::vector<double>>& optional_phase_2_primal_direction,
      const std::optional<ConstraintPartition>& optional_constraint_partition) {
   // form and solve the feasibility problem (with or without constraint partition)
   this->create_current_feasibility_problem(problem, scaling, current_iterate, optional_phase_2_primal_direction, optional_constraint_partition);

   DEBUG << "\nSolving the feasibility subproblem at the current iterate:\n" << current_iterate << "\n";
   Direction feasibility_direction = this->subproblem->solve(statistics, problem, current_iterate);
   feasibility_direction.objective_multiplier = 0.;
   DEBUG << feasibility_direction << "\n";
   assert(feasibility_direction.status == OPTIMAL && "The subproblem was not solved to optimality");

   if (optional_constraint_partition.has_value()) {
      // transfer the constraint partition of the phase-2 direction to the phase-1 direction
      const ConstraintPartition& constraint_partition = optional_constraint_partition.value();
      feasibility_direction.constraint_partition = constraint_partition;
   }
   else {
      // remove the temporary elastic variables
      this->remove_elastic_variables_from_subproblem();
   }
   return feasibility_direction;
}

bool FeasibilityRestoration::is_acceptable(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
      Iterate& trial_iterate, const Direction& direction, PredictedReductionModel& predicted_reduction_model, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed, the progress measures are recomputed\n";
      this->subproblem->subproblem_definition_changed = false;
      this->phase_2_strategy->reset();
      this->subproblem->compute_progress_measures(problem, scaling, current_iterate);
   }

   bool accept = false;
   if (ConstraintRelaxationStrategy::is_small_step(direction)) {
      accept = true;
   }
   else {
      // possibly switch between phase 1 (restoration) and phase 2 (optimality)
      GlobalizationStrategy& current_phase_strategy = this->switch_phase(problem, scaling, current_iterate, trial_iterate, direction);

      // evaluate the predicted reduction
      const double predicted_reduction = predicted_reduction_model.evaluate(step_length);

      // invoke the globalization strategy for acceptance
      accept = current_phase_strategy.check_acceptance(statistics, current_iterate.progress, trial_iterate.progress, direction.objective_multiplier,
            predicted_reduction);
   }

   if (accept) {
      statistics.add_statistic("phase", static_cast<int>(this->current_phase));
      // correct multipliers for infeasibility problem with constraint partition
      if (direction.objective_multiplier == 0. && direction.constraint_partition.has_value()) {
         const ConstraintPartition& constraint_partition = direction.constraint_partition.value();
         FeasibilityRestoration::set_restoration_multipliers(trial_iterate.multipliers.constraints, constraint_partition);
      }
      this->subproblem->compute_residuals(problem, scaling, trial_iterate, direction.objective_multiplier);
   }
   return accept;
}

void FeasibilityRestoration::create_current_feasibility_problem(const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
      const std::optional<std::vector<double>>& optional_phase_2_primal_direction,
      const std::optional<ConstraintPartition>& optional_constraint_partition) {
   // if a constraint partition is given, form a partitioned l1 feasibility problem
   if (optional_constraint_partition.has_value()) {
      const ConstraintPartition& constraint_partition = optional_constraint_partition.value();
      assert(!constraint_partition.infeasible.empty() && "The subproblem is infeasible but no constraint is infeasible");
      // set the multipliers of the violated constraints
      FeasibilityRestoration::set_restoration_multipliers(this->subproblem->constraints_multipliers, constraint_partition);

      // compute the objective model with a zero objective multiplier
      this->subproblem->objective_gradient.clear();
      this->subproblem->build_objective_model(problem, scaling, current_iterate, 0.);

      // assemble the linear objective (sum of the gradients of the violated constraints)
      this->subproblem->compute_feasibility_linear_objective(current_iterate, constraint_partition);

      // update the bounds of the constraints
      this->subproblem->generate_feasibility_bounds(problem, current_iterate.subproblem_constraints, constraint_partition);
   }
   else {
      // no constraint partition given, form an l1 feasibility problem by adding elastic variables
      initialize_vector(this->subproblem->constraints_multipliers, 0.);
      this->subproblem->build_objective_model(problem, scaling, current_iterate, 0.);
      this->add_elastic_variables_to_subproblem(problem, current_iterate);
   }
   // start from the phase-2 solution
   if (optional_phase_2_primal_direction.has_value()) {
      const std::vector<double>& phase_2_primal_direction = optional_phase_2_primal_direction.value();
      this->subproblem->set_initial_point(phase_2_primal_direction);
   }
}

GlobalizationStrategy& FeasibilityRestoration::switch_phase(const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
      Iterate& trial_iterate, const Direction& direction) {
   // possibly go from 1 (restoration) to phase 2 (optimality)
   if (this->current_phase == FEASIBILITY_RESTORATION && 0. < direction.objective_multiplier) {
      // TODO && this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective))
      this->current_phase = OPTIMALITY;
      DEBUG << "Switching from restoration to optimality phase\n";
      if (direction.constraint_partition.has_value()) {
         // remove elastics from the current iterate
         current_iterate.set_number_variables(this->number_subproblem_variables);
      }
      this->evaluate_constraints(problem, scaling, current_iterate);
      this->subproblem->compute_progress_measures(problem, scaling, current_iterate);
   }
   // possibly go from phase 2 (optimality) to 1 (restoration)
   else if (this->current_phase == OPTIMALITY && direction.objective_multiplier == 0.) {
      this->current_phase = FEASIBILITY_RESTORATION;
      DEBUG << "Switching from optimality to restoration phase\n";
      this->phase_2_strategy->notify(current_iterate);
      this->phase_1_strategy->reset();
      this->compute_infeasibility_measures(problem, scaling, current_iterate, direction.constraint_partition);
      this->phase_1_strategy->notify(current_iterate);
   }

   // evaluate the progress measures of the trial iterate
   if (this->current_phase == OPTIMALITY) {
      trial_iterate.set_number_variables(this->number_subproblem_variables);
      this->evaluate_constraints(problem, scaling, trial_iterate);
      this->subproblem->compute_progress_measures(problem, scaling, trial_iterate);
   }
   else { // restoration phase
      this->compute_infeasibility_measures(problem, scaling, trial_iterate, direction.constraint_partition);
   }
   // return the globalization strategy of the current phase
   return (this->current_phase == OPTIMALITY) ? *this->phase_2_strategy : *this->phase_1_strategy;
}

void FeasibilityRestoration::set_restoration_multipliers(std::vector<double>& constraints_multipliers, const ConstraintPartition&
constraint_partition) {
   // the values are derived from the KKT conditions of the feasibility problem
   for (size_t j: constraint_partition.lower_bound_infeasible) {
      constraints_multipliers[j] = 1.;
   }
   for (size_t j: constraint_partition.upper_bound_infeasible) {
      constraints_multipliers[j] = -1.;
   }
   // otherwise, leave the multiplier as it is
}

void FeasibilityRestoration::compute_infeasibility_measures(const Problem& problem, const Scaling& scaling, Iterate& iterate,
      const std::optional<ConstraintPartition>& optional_constraint_partition) {
   // optimality measure: residual of linearly infeasible constraints
   if (optional_constraint_partition.has_value()) {
      const ConstraintPartition& constraint_partition = optional_constraint_partition.value();
      // feasibility measure: residual of all constraints
      this->evaluate_constraints(problem, scaling, iterate);
      const double feasibility_measure = this->subproblem->compute_constraint_violation(problem, scaling, iterate);
      const double objective_measure = this->subproblem->compute_constraint_violation(problem, scaling, iterate, constraint_partition.infeasible);
      iterate.progress = {feasibility_measure, objective_measure};
   }
   else {
      // if no constraint partition is available, simply compute the standard progress measures
      this->evaluate_relaxed_constraints(problem, scaling, iterate);
      this->subproblem->compute_progress_measures(problem, scaling, iterate);
      // add elastic variables to the optimality measure
      this->elastic_variables.negative.for_each_index([&](size_t i) {
         iterate.progress.objective += this->elastic_objective_coefficient*iterate.x[i];
      });
      this->elastic_variables.positive.for_each_index([&](size_t i) {
         iterate.progress.objective += this->elastic_objective_coefficient*iterate.x[i];
      });
   }
}