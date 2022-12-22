// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

FeasibilityRestoration::FeasibilityRestoration(const Model& model, const Options& options) :
      ConstraintRelaxationStrategy(model, options),
      // create the phase-2 optimality problem (= original model)
      optimality_problem(model),
      // create the phase-1 feasibility problem (objective multiplier = 0)
      feasibility_problem(model, 0.),
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

   // compute the progress measures of the initial point
   this->set_infeasibility_measure(first_iterate);
   this->set_scaled_optimality_measure(first_iterate);
   this->subproblem->set_unscaled_optimality_measure(this->get_current_reformulated_problem(), first_iterate);

   // compute the residuals of the initial point
   ConstraintRelaxationStrategy::evaluate_reformulation_functions(this->optimality_problem, first_iterate);
   this->compute_optimality_condition_residuals(this->optimality_problem, first_iterate);

   // initialize the globalization strategies
   this->phase_1_strategy->initialize(first_iterate);
   this->phase_2_strategy->initialize(first_iterate);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) {
   DEBUG << "Current iterate\n" << current_iterate << '\n';
   if (this->current_phase == Phase::OPTIMALITY) {
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
   if (direction.status == SubproblemStatus::INFEASIBLE) {
      direction = this->solve_feasibility_problem(statistics, current_iterate, direction.primals);
   }
   return direction;
}

// form and solve the feasibility problem (with an initial point)
Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate, const std::vector<double>& initial_point) {
   this->subproblem->set_initial_point(initial_point);
   return this->solve_feasibility_problem(statistics, current_iterate);
}

// form and solve the feasibility problem
Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate) {
   this->subproblem->prepare_for_feasibility_problem(current_iterate);

   // set the initial values of the elastic variables
   this->subproblem->set_elastic_variable_values(this->feasibility_problem, current_iterate);

   DEBUG << "Solving the feasibility subproblem\n";
   Direction direction = this->subproblem->solve(statistics, this->feasibility_problem, current_iterate);
   direction.objective_multiplier = 0.;
   direction.norm = norm_inf(direction.primals, Range(this->optimality_problem.number_variables));
   DEBUG << direction << '\n';
   assert(direction.status == SubproblemStatus::OPTIMAL && "The feasibility subproblem was not solved to optimality");
   return direction;
}

void FeasibilityRestoration::compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction) {
   // refresh the unscaled optimality measures for the current iterate
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed, the optimality measure is recomputed\n";
      this->subproblem->set_unscaled_optimality_measure(this->get_current_reformulated_problem(), current_iterate);
      this->phase_2_strategy->reset();
      this->subproblem->subproblem_definition_changed = false;
   }

   // possibly go from optimality phase to restoration phase
   if (this->current_phase == Phase::OPTIMALITY && direction.objective_multiplier == 0.) {
      this->switch_to_feasibility_restoration(current_iterate);
   }
   // possibly go from restoration phase to optimality phase
   else if (this->current_phase == Phase::FEASIBILITY_RESTORATION) {
      // evaluate measure of infeasibility ("scaled optimality" quantity in phase-1 definition)
      this->set_scaled_optimality_measure(trial_iterate);
      if (ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model, current_iterate, direction, 1.) == 0. &&
            this->phase_2_strategy->is_feasibility_iterate_acceptable(trial_iterate.nonlinear_progress.scaled_optimality)) {
         this->switch_to_optimality(current_iterate, trial_iterate);
      }
   }

   // evaluate the progress measures of the trial iterate
   this->set_infeasibility_measure(trial_iterate);
   this->set_scaled_optimality_measure(trial_iterate);
   this->subproblem->set_unscaled_optimality_measure(this->get_current_reformulated_problem(), trial_iterate);
}

void FeasibilityRestoration::switch_to_feasibility_restoration(Iterate& current_iterate) {
   DEBUG << "Switching from optimality to restoration phase\n";
   this->current_phase = Phase::FEASIBILITY_RESTORATION;
   this->phase_2_strategy->register_current_progress(current_iterate.nonlinear_progress);
   this->phase_1_strategy->reset();
   // refresh the progress measures of the current iterate
   this->set_scaled_optimality_measure(current_iterate);
   this->set_infeasibility_measure(current_iterate);
   this->phase_1_strategy->register_current_progress(current_iterate.nonlinear_progress);
}

void FeasibilityRestoration::switch_to_optimality(Iterate& current_iterate, Iterate& trial_iterate) {
   DEBUG << "Switching from restoration to optimality phase\n";
   this->current_phase = Phase::OPTIMALITY;
   current_iterate.set_number_variables(this->optimality_problem.number_variables);
   trial_iterate.set_number_variables(this->optimality_problem.number_variables);
   // refresh the progress measures of current iterate
   this->set_scaled_optimality_measure(current_iterate);
   this->set_infeasibility_measure(current_iterate);
}

bool FeasibilityRestoration::is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      double step_length) {
   this->compute_progress_measures(current_iterate, trial_iterate, direction);

   bool accept = false;
   if (this->is_small_step(direction)) {
      DEBUG << "Small step acceptable\n";
      accept = true;
   }
   else {
      // evaluate the predicted reduction
      PredictedReductionModels predicted_reduction_model = this->generate_predicted_reduction_models(current_iterate, direction);
      DEBUG << "Infeasibility model:       " << predicted_reduction_model.infeasibility.text << '\n';
      DEBUG << "Scaled optimality model:   " << predicted_reduction_model.scaled_optimality.text << '\n';
      DEBUG << "Unscaled optimality model: " << predicted_reduction_model.unscaled_optimality.text << '\n';
      ProgressMeasures predicted_reduction = {
            predicted_reduction_model.infeasibility(step_length),
            predicted_reduction_model.scaled_optimality(step_length),
            predicted_reduction_model.unscaled_optimality(step_length)
      };
      // invoke the globalization strategy for acceptance
      GlobalizationStrategy& current_phase_strategy = this->get_current_globalization_strategy();
      accept = current_phase_strategy.is_iterate_acceptable(current_iterate.nonlinear_progress, trial_iterate.nonlinear_progress, predicted_reduction);
   }
   if (accept) {
      statistics.add_statistic("phase", static_cast<int>(this->current_phase));
      ConstraintRelaxationStrategy::evaluate_reformulation_functions(this->get_current_reformulated_problem(), trial_iterate);
      this->compute_optimality_condition_residuals(this->get_current_reformulated_problem(), trial_iterate);
   }
   return accept;
}

PredictedReductionModels FeasibilityRestoration::generate_predicted_reduction_models(const Iterate& current_iterate, const Direction& direction) const {
   return {
      this->generate_predicted_infeasibility_reduction_model(current_iterate, direction),
      this->generate_predicted_scaled_optimality_reduction_model(current_iterate, direction),
      this->subproblem->generate_predicted_unscaled_optimality_reduction_model(this->get_current_reformulated_problem(), current_iterate, direction)
   };
}

const NonlinearProblem& FeasibilityRestoration::get_current_reformulated_problem() const {
   if (this->current_phase == Phase::OPTIMALITY) {
      return this->optimality_problem;
   }
   else {
      return this->feasibility_problem;
   }
}

GlobalizationStrategy& FeasibilityRestoration::get_current_globalization_strategy() const {
   return (this->current_phase == Phase::OPTIMALITY) ? *this->phase_2_strategy : *this->phase_1_strategy;
}

void FeasibilityRestoration::set_variable_bounds(const Iterate& current_iterate, double trust_region_radius) {
   this->subproblem->set_variable_bounds(this->feasibility_problem, current_iterate, trust_region_radius);
}

Direction FeasibilityRestoration::compute_second_order_correction(Iterate& trial_iterate) {
   return this->subproblem->compute_second_order_correction(this->get_current_reformulated_problem(), trial_iterate);
}

void FeasibilityRestoration::set_infeasibility_measure(Iterate& iterate) {
   if (this->current_phase == Phase::OPTIMALITY) {
      // constraint violation
      iterate.evaluate_constraints(this->original_model);
      iterate.nonlinear_progress.infeasibility = this->original_model.compute_constraint_violation(iterate.model_evaluations.constraints, L1_NORM);
   }
   else {
      // 0
      iterate.nonlinear_progress.infeasibility = 0.;
   }
}

PredictedReductionModel FeasibilityRestoration::generate_predicted_infeasibility_reduction_model(const Iterate& current_iterate,
      const Direction& direction) const {
   if (this->current_phase == Phase::OPTIMALITY) {
      return {[&](double step_length) {
         const double linearized_constraint_violation = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model,
               current_iterate, direction, step_length);
         return current_iterate.primal_constraint_violation - linearized_constraint_violation;
      }, "‖c(x)‖₁ - ‖c(x) + ∇c(x)^T (αd)‖₁"};
   }
   else {
      return {[](double /*step_length*/) {
         return 0.;
      }, "0"};
   }
}

void FeasibilityRestoration::set_scaled_optimality_measure(Iterate& iterate) {
   if (this->current_phase == Phase::OPTIMALITY) {
      // original objective
      iterate.evaluate_objective(this->original_model);
      iterate.nonlinear_progress.scaled_optimality = iterate.model_evaluations.objective;
   }
   else {
      // constraint violation
      iterate.evaluate_constraints(this->original_model);
      iterate.nonlinear_progress.scaled_optimality = this->original_model.compute_constraint_violation(iterate.model_evaluations.constraints, L1_NORM);
   }
}

PredictedReductionModel FeasibilityRestoration::generate_predicted_scaled_optimality_reduction_model(const Iterate& current_iterate,
      const Direction& direction) const {
   if (this->current_phase == Phase::OPTIMALITY) {
      // precompute expensive quantities
      const double scaled_directional_derivative = dot(direction.primals, current_iterate.model_evaluations.objective_gradient);
      return {[=](double step_length) {
         // return a function of the step length that cheaply assembles the predicted reduction
         return step_length * (-scaled_directional_derivative);
      }, "-∇f(x)^T (αd)"};
   }
   else {
      return {[&](double step_length) {
         const double linearized_constraint_violation = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model,
               current_iterate, direction, step_length);
         return current_iterate.primal_constraint_violation - linearized_constraint_violation;
      }, "‖c(x)‖₁ - ‖c(x) + ∇c(x)^T (αd)‖₁"};
   }
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
