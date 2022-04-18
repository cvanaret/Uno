#include <cassert>
#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

FeasibilityRestoration::FeasibilityRestoration(const Model& model, const Options& options) :
      // create the phase-1 feasibility problem (objective multiplier = 0) with elastic variables
      ConstraintRelaxationStrategy(model, 0., options),
      // create the globalization strategies (one for each phase)
      phase_1_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)),
      phase_2_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)) {
}

void FeasibilityRestoration::initialize(Statistics& statistics, Iterate& first_iterate) {
   statistics.add_column("phase", Statistics::int_width, 4);

   // initialize the subproblem
   this->subproblem->initialize(statistics, this->optimality_problem, first_iterate);

   // compute the progress measures and the residuals of the initial point
   first_iterate.nonlinear_progress.infeasibility = this->compute_infeasibility_measure(first_iterate);
   first_iterate.nonlinear_progress.objective = this->subproblem->compute_optimality_measure(this->optimality_problem.model, first_iterate);
   this->subproblem->compute_nonlinear_residuals(this->optimality_problem, first_iterate);

   // initialize the globalization strategies
   this->phase_1_strategy->initialize(statistics, first_iterate);
   this->phase_2_strategy->initialize(statistics, first_iterate);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) {
   DEBUG << "Current iterate\n" << current_iterate << "\n";
   if (this->current_phase == OPTIMALITY) {
      return this->solve_optimality_problem(statistics, current_iterate);
   }
   else {
      return this->solve_feasibility_problem(statistics, current_iterate, std::nullopt);
   }
}

Direction FeasibilityRestoration::solve_optimality_problem(Statistics& statistics, Iterate& current_iterate) {
   // solve the subproblem
   DEBUG << "Solving the optimality subproblem\n";
   Direction direction = this->subproblem->solve(statistics, this->optimality_problem, current_iterate);
   direction.objective_multiplier = 1.;
   direction.norm = norm_inf(direction.x, 0, this->optimality_problem.number_variables);
   DEBUG << direction << "\n";

   // infeasible subproblem: try to minimize the constraint violation by solving the feasibility subproblem
   if (direction.status == INFEASIBLE) {
      direction = this->solve_feasibility_problem(statistics, current_iterate, direction.x);
   }
   return direction;
}

// form and solve the feasibility problem (with or without constraint partition)
Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
      const std::optional<std::vector<double>>& optional_phase_2_solution) {
   // register the proximal coefficient and reference point
   this->relaxed_problem.set_proximal_coefficient(this->subproblem->get_proximal_coefficient());
   this->relaxed_problem.set_proximal_reference_point(current_iterate.x);

   // build the objective model of the feasibility problem
   this->subproblem->set_elastic_variables(this->relaxed_problem, current_iterate);

   // start from the phase-2 solution
   this->subproblem->set_initial_point(optional_phase_2_solution);

   DEBUG << "Solving the feasibility subproblem\n";
   Direction feasibility_direction = this->subproblem->solve(statistics, this->relaxed_problem, current_iterate);
   feasibility_direction.objective_multiplier = 0.;
   feasibility_direction.norm = norm_inf(feasibility_direction.x, 0, this->optimality_problem.number_variables);
   // create constraint partition
   ConstraintPartition constraint_partition(this->optimality_problem.number_constraints);
   constraint_partition.infeasible = this->relaxed_problem.get_violated_linearized_constraints(feasibility_direction.x);
   feasibility_direction.constraint_partition = constraint_partition;
   DEBUG << feasibility_direction << "\n";
   assert(feasibility_direction.status == OPTIMAL && "The subproblem was not solved to optimality");
   return feasibility_direction;
}

bool FeasibilityRestoration::is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      PredictedReductionModel& predicted_reduction_model, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed, the optimality measure is recomputed\n";
      current_iterate.nonlinear_progress.objective = this->subproblem->compute_optimality_measure(this->optimality_problem.model, current_iterate);
      this->phase_2_strategy->reset();
      this->subproblem->subproblem_definition_changed = false;
   }

   // possibly switch between phase 1 (restoration) and phase 2 (optimality)
   GlobalizationStrategy& current_phase_strategy = this->switch_phase(current_iterate, trial_iterate, direction);

   bool accept = false;
   if (ConstraintRelaxationStrategy::is_small_step(direction)) {
      accept = true;
   }
   else {
      // evaluate the predicted reduction
      const double predicted_reduction = predicted_reduction_model.evaluate(step_length);

      // invoke the globalization strategy for acceptance
      accept = current_phase_strategy.is_acceptable(statistics, current_iterate.nonlinear_progress, trial_iterate.nonlinear_progress,
            direction.objective_multiplier, predicted_reduction);
   }

   if (accept) {
      statistics.add_statistic("phase", static_cast<int>(this->current_phase));
      if (this->current_phase == OPTIMALITY) {
         this->subproblem->compute_nonlinear_residuals(this->optimality_problem, trial_iterate);
      }
      else {
         this->subproblem->compute_nonlinear_residuals(this->relaxed_problem, trial_iterate);
      }
      trial_iterate.evaluate_objective(this->optimality_problem.model);
   }
   return accept;
}

GlobalizationStrategy& FeasibilityRestoration::switch_phase(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction) {
   if (this->current_phase == OPTIMALITY && direction.objective_multiplier == 0.) {
      this->current_phase = FEASIBILITY_RESTORATION;
      DEBUG << "Switching from optimality to restoration phase\n";
      this->phase_2_strategy->notify(current_iterate);
      this->phase_1_strategy->reset();
      current_iterate.reset_evaluations();
      // update the measure of optimality
      const std::vector<size_t>& infeasible_constraints = direction.constraint_partition.value().infeasible;
      current_iterate.nonlinear_progress.objective = this->compute_optimality_measure(current_iterate, infeasible_constraints);
      this->phase_1_strategy->notify(current_iterate);
   }
   // possibly go from 1 (restoration) to phase 2 (optimality)
   else if (this->current_phase == FEASIBILITY_RESTORATION && direction.constraint_partition.value().infeasible.empty()) {
      // TODO && this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective))
      this->current_phase = OPTIMALITY;
      DEBUG << "Switching from restoration to optimality phase\n";
      current_iterate.set_number_variables(this->optimality_problem.number_variables);
      current_iterate.reset_evaluations();
      current_iterate.nonlinear_progress.objective = this->subproblem->compute_optimality_measure(this->optimality_problem.model, current_iterate);
      trial_iterate.set_number_variables(this->optimality_problem.number_variables);
   }

   // evaluate the progress measures of the trial iterate
   trial_iterate.nonlinear_progress.infeasibility = this->compute_infeasibility_measure(trial_iterate);
   if (this->current_phase == OPTIMALITY) {
      trial_iterate.nonlinear_progress.objective = this->subproblem->compute_optimality_measure(this->optimality_problem.model, trial_iterate);
   }
   else {
      trial_iterate.nonlinear_progress.objective = this->compute_optimality_measure(trial_iterate, direction.constraint_partition.value().infeasible);
   }

   // return the globalization strategy of the current phase
   return (this->current_phase == OPTIMALITY) ? *this->phase_2_strategy : *this->phase_1_strategy;
}

void FeasibilityRestoration::set_variable_bounds(const Iterate& current_iterate, double trust_region_radius) {
   // set the bounds of all the variables (primal + elastics)
   this->subproblem->set_variable_bounds(this->relaxed_problem.model, current_iterate, trust_region_radius);
}

Direction FeasibilityRestoration::compute_second_order_correction(Iterate& trial_iterate) {
   if (this->current_phase == OPTIMALITY) {
      return this->subproblem->compute_second_order_correction(this->optimality_problem, trial_iterate);
   }
   else {
      return this->subproblem->compute_second_order_correction(this->relaxed_problem, trial_iterate);
   }
}

PredictedReductionModel FeasibilityRestoration::generate_predicted_reduction_model(const Direction& direction) const {
   if (this->current_phase == OPTIMALITY) {
      return this->subproblem->generate_predicted_reduction_model(this->optimality_problem, direction);
   }
   else {
      return this->subproblem->generate_predicted_reduction_model(this->relaxed_problem, direction);
   }
}

double FeasibilityRestoration::compute_infeasibility_measure(Iterate& iterate) {
   iterate.evaluate_constraints(this->optimality_problem.model);
   return this->optimality_problem.compute_constraint_violation(iterate.original_evaluations.constraints, L1_NORM);
}

double FeasibilityRestoration::compute_optimality_measure(Iterate& iterate, const std::vector<size_t>& infeasible_constraints) {
   iterate.evaluate_constraints(this->optimality_problem.model);
   return this->optimality_problem.compute_constraint_violation(iterate.original_evaluations.constraints, infeasible_constraints, L1_NORM);
}

void FeasibilityRestoration::register_accepted_iterate(Iterate& iterate) {
   if (this->current_phase == OPTIMALITY) {
      this->subproblem->postprocess_accepted_iterate(this->optimality_problem, iterate);
   }
   else {
      this->subproblem->postprocess_accepted_iterate(this->relaxed_problem, iterate);
   }
}