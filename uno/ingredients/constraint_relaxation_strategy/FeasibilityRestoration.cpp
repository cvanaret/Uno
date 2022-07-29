// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

FeasibilityRestoration::FeasibilityRestoration(const Model& model, const Options& options) :
      ConstraintRelaxationStrategy(model, options),
      // create the optimality problem
      optimality_problem(model),
      // create the phase-1 feasibility problem (objective multiplier = 0)
      feasibility_problem(model, 0., options.get_double("l1_constraint_violation_coefficient"), options.get_bool("l1_use_proximal_term")),
      subproblem(SubproblemFactory::create(this->feasibility_problem.number_variables, this->feasibility_problem.number_constraints,
            this->feasibility_problem.get_maximum_number_hessian_nonzeros(), options)),
      // create the globalization strategies (one for each phase)
      phase_1_strategy(GlobalizationStrategyFactory::create(options.get_string("strategy"), options)),
      phase_2_strategy(GlobalizationStrategyFactory::create(options.get_string("strategy"), options)),
      statistics_restoration_phase_column_order(options.get_int("statistics_restoration_phase_column_order")) {
}

void FeasibilityRestoration::initialize(Statistics& statistics, Iterate& first_iterate) {
   statistics.add_column("phase", Statistics::int_width, this->statistics_restoration_phase_column_order);

   // initialize the subproblem
   this->subproblem->initialize(statistics, this->optimality_problem, first_iterate);

   // compute the progress measures and the residuals of the initial point
   first_iterate.nonlinear_progress.infeasibility = this->compute_infeasibility_measure(first_iterate);
   first_iterate.nonlinear_progress.optimality = this->subproblem->compute_optimality_measure(this->optimality_problem, first_iterate);
   this->compute_nonlinear_residuals(this->optimality_problem, first_iterate);

   // initialize the globalization strategies
   this->phase_1_strategy->initialize(first_iterate);
   this->phase_2_strategy->initialize(first_iterate);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) {
   DEBUG << "Current iterate\n" << current_iterate << '\n';
   if (this->current_phase == OPTIMALITY) {
      return this->solve_optimality_problem(statistics, current_iterate);
   }
   else {
      return this->solve_feasibility_problem(statistics, current_iterate);
   }
}

Direction FeasibilityRestoration::solve_optimality_problem(Statistics& statistics, Iterate& current_iterate) {
   // solve the subproblem
   DEBUG << "Solving the optimality subproblem\n";
   Direction direction = this->subproblem->solve(statistics, this->optimality_problem, current_iterate);
   direction.objective_multiplier = 1.;
   direction.norm = norm_inf(direction.primals, Range(this->optimality_problem.number_variables));
   DEBUG << direction << '\n';

   // infeasible subproblem: try to minimize the constraint violation by solving the feasibility subproblem
   if (direction.status == Status::INFEASIBLE) {
      direction = this->solve_feasibility_problem(statistics, current_iterate, direction.primals);
   }
   return direction;
}

// form and solve the feasibility problem (with or without constraint partition)
Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate) {
   // register the proximal coefficient and reference point
   this->feasibility_problem.set_proximal_coefficient(this->subproblem->get_proximal_coefficient());
   this->feasibility_problem.set_proximal_reference_point(current_iterate.primals);

   // build the objective model of the feasibility problem
   this->subproblem->set_elastic_variables(this->feasibility_problem, current_iterate);

   DEBUG << "Solving the feasibility subproblem\n";
   Direction feasibility_direction = this->subproblem->solve(statistics, this->feasibility_problem, current_iterate);
   assert(feasibility_direction.status == Status::OPTIMAL && "The feasibility subproblem was not solved to optimality");
   feasibility_direction.objective_multiplier = 0.;
   feasibility_direction.norm = norm_inf(feasibility_direction.primals, Range(this->optimality_problem.number_variables));
   DEBUG << feasibility_direction << '\n';
   return feasibility_direction;
}

// form and solve the feasibility problem (with or without constraint partition)
Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
      const std::vector<double>& initial_point) {
   this->subproblem->set_initial_point(initial_point);
   return this->solve_feasibility_problem(statistics, current_iterate);
}

bool FeasibilityRestoration::is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      PredictedOptimalityReductionModel& predicted_optimality_reduction_model, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed, the optimality measure is recomputed\n";
      current_iterate.nonlinear_progress.optimality = this->subproblem->compute_optimality_measure(this->optimality_problem, current_iterate);
      this->phase_2_strategy->reset();
      this->subproblem->subproblem_definition_changed = false;
   }

   // possibly switch between phase 1 (restoration) and phase 2 (optimality)
   GlobalizationStrategy& current_phase_strategy = this->switch_phase(current_iterate, trial_iterate, direction);

   bool accept = false;
   if (this->is_small_step(direction)) {
      accept = true;
   }
   else {
      // evaluate the predicted reduction
      const ProgressMeasures predicted_reduction = {
            ConstraintRelaxationStrategy::compute_predicted_infeasibility_reduction(this->model, current_iterate, direction, step_length),
            predicted_optimality_reduction_model.evaluate(step_length)
      };

      // invoke the globalization strategy for acceptance
      accept = current_phase_strategy.is_acceptable(current_iterate.nonlinear_progress, trial_iterate.nonlinear_progress,
            direction.objective_multiplier, predicted_reduction);
   }

   if (accept) {
      statistics.add_statistic("phase", static_cast<int>(this->current_phase));
      this->compute_nonlinear_residuals(this->get_current_reformulated_problem(), trial_iterate);
   }
   return accept;
}

PredictedOptimalityReductionModel FeasibilityRestoration::generate_predicted_optimality_reduction_model(const Direction& direction) const {
   return this->subproblem->generate_predicted_optimality_reduction_model(this->get_current_reformulated_problem(), direction);
}

const NonlinearProblem& FeasibilityRestoration::get_current_reformulated_problem() const {
   if (this->current_phase == OPTIMALITY) {
      return this->optimality_problem;
   }
   else {
      return this->feasibility_problem;
   }
}

GlobalizationStrategy& FeasibilityRestoration::switch_phase(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction) {
   std::vector<size_t> infeasible_constraints = this->feasibility_problem.get_violated_linearized_constraints(direction.primals);

   if (this->current_phase == OPTIMALITY && direction.objective_multiplier == 0.) {
      this->switch_to_feasibility_restoration(current_iterate, infeasible_constraints);
   }
   // possibly go from 1 (restoration) to phase 2 (optimality)
   else if (this->current_phase == FEASIBILITY_RESTORATION && infeasible_constraints.empty()) {
      // TODO && this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective))
      this->switch_to_optimality(current_iterate, trial_iterate);
   }

   // evaluate the progress measures of the trial iterate
   trial_iterate.evaluate_objective(this->model);
   trial_iterate.nonlinear_progress.infeasibility = this->compute_infeasibility_measure(trial_iterate);
   if (this->current_phase == OPTIMALITY) {
      trial_iterate.nonlinear_progress.optimality = this->subproblem->compute_optimality_measure(this->optimality_problem, trial_iterate);
   }
   else {
      // TODO check that
      current_iterate.nonlinear_progress.optimality = this->compute_optimality_measure(current_iterate, infeasible_constraints);
      trial_iterate.nonlinear_progress.optimality = this->compute_optimality_measure(trial_iterate, infeasible_constraints);
   }

   // return the globalization strategy of the current phase
   return (this->current_phase == OPTIMALITY) ? *this->phase_2_strategy : *this->phase_1_strategy;
}

void FeasibilityRestoration::switch_to_feasibility_restoration(Iterate& current_iterate, const std::vector<size_t>& infeasible_constraints) {
   this->current_phase = FEASIBILITY_RESTORATION;
   DEBUG << "Switching from optimality to restoration phase\n";
   this->phase_2_strategy->notify(current_iterate);
   this->phase_1_strategy->reset();
   // update the measure of optimality
   current_iterate.nonlinear_progress.optimality = this->compute_optimality_measure(current_iterate, infeasible_constraints);
   this->phase_1_strategy->notify(current_iterate);
}

void FeasibilityRestoration::switch_to_optimality(Iterate& current_iterate, Iterate& trial_iterate) {
   this->current_phase = OPTIMALITY;
   DEBUG << "Switching from restoration to optimality phase\n";
   current_iterate.set_number_variables(this->optimality_problem.number_variables);
   current_iterate.nonlinear_progress.optimality = this->subproblem->compute_optimality_measure(this->optimality_problem, current_iterate);
   trial_iterate.set_number_variables(this->optimality_problem.number_variables);
}

void FeasibilityRestoration::set_variable_bounds(const Iterate& current_iterate, double trust_region_radius) {
   // set the bounds of all the variables (primal + elastics)
   this->subproblem->set_variable_bounds(this->feasibility_problem, current_iterate, trust_region_radius);
}

Direction FeasibilityRestoration::compute_second_order_correction(Iterate& trial_iterate) {
   return this->subproblem->compute_second_order_correction(this->get_current_reformulated_problem(), trial_iterate);
}

double FeasibilityRestoration::compute_infeasibility_measure(Iterate& iterate) {
   iterate.evaluate_constraints(this->model);
   return this->model.compute_constraint_violation(iterate.original_evaluations.constraints, L1_NORM);
}

double FeasibilityRestoration::compute_optimality_measure(Iterate& iterate, const std::vector<size_t>& infeasible_constraints) {
   iterate.evaluate_constraints(this->model);
   return this->model.compute_constraint_violation(iterate.original_evaluations.constraints, infeasible_constraints, L1_NORM);
}

void FeasibilityRestoration::register_accepted_iterate(Iterate& iterate) {
   this->subproblem->postprocess_accepted_iterate(this->get_current_reformulated_problem(), iterate);
}

size_t FeasibilityRestoration::get_hessian_evaluation_count() const {
   return this->subproblem->get_hessian_evaluation_count();
}

size_t FeasibilityRestoration::get_number_subproblems_solved() const {
   return this->subproblem->number_subproblems_solved;
}