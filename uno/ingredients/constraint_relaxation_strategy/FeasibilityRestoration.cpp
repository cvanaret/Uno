// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"

FeasibilityRestoration::FeasibilityRestoration(Statistics& statistics, const Model& model, const Options& options) :
      ConstraintRelaxationStrategy(model, options),
      // create the (optimality phase) optimality problem (= original model)
      optimality_problem(model),
      // create the (restoration phase) feasibility problem (objective multiplier = 0)
      feasibility_problem(model, 0., options.get_double("l1_constraint_violation_coefficient")),
      subproblem(SubproblemFactory::create(statistics, this->feasibility_problem.number_variables, this->feasibility_problem.number_constraints,
            this->feasibility_problem.get_number_jacobian_nonzeros(), this->feasibility_problem.get_number_hessian_nonzeros(), options)),
      // create the globalization strategies (one for each phase)
      restoration_phase_strategy(GlobalizationStrategyFactory::create(statistics, options.get_string("globalization_strategy"), false, options)),
      optimality_phase_strategy(GlobalizationStrategyFactory::create(statistics, options.get_string("globalization_strategy"), true, options)),
      l1_constraint_violation_coefficient(options.get_double("l1_constraint_violation_coefficient")),
      tolerance(options.get_double("tolerance")),
      test_linearized_feasibility(options.get_bool("feasibility_restoration_test_linearized_feasibility")) {
   statistics.add_column("phase", Statistics::int_width, options.get_int("statistics_restoration_phase_column_order"));
}

void FeasibilityRestoration::initialize(Iterate& initial_iterate) {
   this->subproblem->generate_initial_iterate(this->optimality_problem, initial_iterate);

   // compute the progress measures and residuals of the initial point
   this->set_progress_measures_for_optimality_problem(initial_iterate);
   ConstraintRelaxationStrategy::compute_primal_dual_residuals(this->optimality_problem, initial_iterate, this->residual_norm);

   // initialize the globalization strategies
   this->restoration_phase_strategy->initialize(initial_iterate);
   this->optimality_phase_strategy->initialize(initial_iterate);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
      WarmstartInformation& warmstart_information) {
   // solve the optimality problem
   if (this->current_phase == Phase::OPTIMALITY) {
      try {
         DEBUG << "Solving the optimality subproblem\n";
         Direction direction = this->solve_subproblem(statistics, this->optimality_problem, current_iterate, warmstart_information);
         // infeasible subproblem: switch to the feasibility problem, starting from the current direction
         if (direction.status == SubproblemStatus::INFEASIBLE) {
            this->switch_to_feasibility_problem(current_iterate, warmstart_information);
            this->subproblem->set_initial_point(direction.primals);
         }
         else {
            // things ran smoothly: return the direction
            return direction;
         }
      }
      catch (const UnstableRegularization&) {
         this->switch_to_feasibility_problem(current_iterate, warmstart_information);
      }
   }

   // feasibility problem: minimize constraint violation
   DEBUG << "Solving the feasibility subproblem\n";
   // note: failure of regularization should not happen here, since the feasibility Jacobian is full rank
   return this->solve_subproblem(statistics, this->feasibility_problem, current_iterate, warmstart_information);
}

// an initial point is provided
Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
      const std::vector<double>& initial_point, WarmstartInformation& warmstart_information) {
   this->subproblem->set_initial_point(initial_point);
   return this->compute_feasible_direction(statistics, current_iterate, warmstart_information);
}

void FeasibilityRestoration::switch_to_feasibility_problem(Iterate& current_iterate, WarmstartInformation& warmstart_information) {
   if (this->current_phase == Phase::FEASIBILITY_RESTORATION) {
      throw std::runtime_error("NOOOOPE\n");
   }
   this->switch_to_feasibility_restoration(current_iterate, warmstart_information);
}

Direction FeasibilityRestoration::solve_subproblem(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
      WarmstartInformation& warmstart_information) {
   if (this->switched_to_optimality_phase) {
      this->switched_to_optimality_phase = false;
      warmstart_information.set_cold_start();
   }

   Direction direction = this->subproblem->solve(statistics, problem, current_iterate, warmstart_information);
   direction.norm = norm_inf(direction.primals, Range(this->optimality_problem.number_variables));
   direction.multipliers.objective = problem.get_objective_multiplier();
   DEBUG2 << direction << '\n';
   return direction;
}

void FeasibilityRestoration::compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      double step_length) {
   // refresh the auxiliary measure for the current iterate
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed, the auxiliary measure is recomputed\n";
      this->restoration_phase_strategy->reset();
      this->optimality_phase_strategy->reset();
      this->subproblem->set_auxiliary_measure(this->current_problem(), current_iterate);
      this->subproblem->subproblem_definition_changed = false;
   }

   // possibly go from restoration phase to optimality phase
   if (this->current_phase == Phase::FEASIBILITY_RESTORATION && (not this->test_linearized_feasibility ||
         ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model, current_iterate, direction, step_length) <=
               this->tolerance)) {
      // if the trial infeasibility improves upon the best known infeasibility of the globalization strategy
      trial_iterate.evaluate_constraints(this->original_model);
      const double trial_infeasibility = this->original_model.compute_constraint_violation(trial_iterate.evaluations.constraints,
            this->progress_norm);
      if (this->optimality_phase_strategy->is_infeasibility_acceptable(trial_infeasibility)) {
         this->switch_to_optimality(current_iterate, trial_iterate);
      }
   }

   // evaluate the progress measures of the trial iterate
   if (this->current_phase == Phase::OPTIMALITY) {
      this->set_progress_measures_for_optimality_problem(trial_iterate);
   }
   else {
      this->set_progress_measures_for_feasibility_problem(trial_iterate);
   }
}

void FeasibilityRestoration::switch_to_feasibility_restoration(Iterate& current_iterate, WarmstartInformation& warmstart_information) {
   DEBUG << "Switching from optimality to restoration phase\n";
   this->current_phase = Phase::FEASIBILITY_RESTORATION;
   this->optimality_phase_strategy->register_current_progress(current_iterate.progress);
   this->subproblem->initialize_feasibility_problem();
   this->subproblem->set_elastic_variable_values(this->feasibility_problem, current_iterate);
   DEBUG2 << "Current iterate:\n" << current_iterate << '\n';

   // refresh the progress measures of the current iterate
   this->set_progress_measures_for_feasibility_problem(current_iterate);

   current_iterate.multipliers.objective = 0.;
   this->restoration_phase_strategy->reset();
   this->restoration_phase_strategy->register_current_progress(current_iterate.progress);
   warmstart_information.set_cold_start();
}

void FeasibilityRestoration::switch_to_optimality(Iterate& current_iterate, Iterate& trial_iterate) {
   DEBUG << "Switching from restoration to optimality phase\n";
   this->current_phase = Phase::OPTIMALITY;
   current_iterate.set_number_variables(this->optimality_problem.number_variables);
   trial_iterate.set_number_variables(this->optimality_problem.number_variables);
   this->subproblem->exit_feasibility_problem(this->optimality_problem, trial_iterate);
   this->switched_to_optimality_phase = true;

   // refresh the progress measures of current iterate
   this->set_progress_measures_for_optimality_problem(current_iterate);
   current_iterate.multipliers.objective = 1.;
   trial_iterate.multipliers.objective = 1.;
}

bool FeasibilityRestoration::is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      double step_length) {
   this->subproblem->postprocess_iterate(this->current_problem(), trial_iterate);
   this->compute_progress_measures(current_iterate, trial_iterate, direction, step_length);

   bool accept_iterate = false;
   if (direction.norm == 0.) {
      DEBUG << "Zero step acceptable\n";
      trial_iterate.evaluate_objective(this->original_model);
      accept_iterate = true;
   }
   else {
      // evaluate the predicted reduction
      ProgressMeasures predicted_reduction = (this->current_phase == Phase::OPTIMALITY) ?
            this->compute_predicted_reduction_models_for_optimality_problem(current_iterate, direction, step_length) :
            this->compute_predicted_reduction_models_for_feasibility_problem(current_iterate, direction, step_length);

      // invoke the globalization strategy for acceptance
      GlobalizationStrategy& current_phase_strategy = this->current_globalization_strategy();
      accept_iterate = current_phase_strategy.is_iterate_acceptable(statistics, trial_iterate, current_iterate.progress, trial_iterate.progress,
            predicted_reduction, this->current_problem().get_objective_multiplier());
   }

   if (accept_iterate) {
      // compute the primal-dual residuals
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(this->optimality_problem, trial_iterate, this->residual_norm);
      if (this->current_phase == Phase::OPTIMALITY) {
         statistics.add_statistic("complementarity", trial_iterate.residuals.optimality_complementarity);
         statistics.add_statistic("stationarity", trial_iterate.residuals.optimality_stationarity);
         if (this->original_model.is_constrained()) {
            statistics.add_statistic("primal infeas.", trial_iterate.progress.infeasibility);
         }
      }
      else {
         statistics.add_statistic("complementarity", trial_iterate.residuals.feasibility_complementarity);
         statistics.add_statistic("stationarity", trial_iterate.residuals.feasibility_stationarity);
         if (this->original_model.is_constrained()) {
            statistics.add_statistic("primal infeas.", trial_iterate.progress.optimality(1.));
         }
      }
      statistics.add_statistic("phase", static_cast<int>(this->current_phase));
   }
   return accept_iterate;
}

const NonlinearProblem& FeasibilityRestoration::current_problem() const {
   if (this->current_phase == Phase::OPTIMALITY) {
      return this->optimality_problem;
   }
   else {
      return this->feasibility_problem;
   }
}

GlobalizationStrategy& FeasibilityRestoration::current_globalization_strategy() const {
   return (this->current_phase == Phase::OPTIMALITY) ? *this->optimality_phase_strategy : *this->restoration_phase_strategy;
}

void FeasibilityRestoration::set_trust_region_radius(double trust_region_radius) {
   this->subproblem->set_trust_region_radius(trust_region_radius);
}

/* progress measures */

void FeasibilityRestoration::set_progress_measures_for_optimality_problem(Iterate& iterate) {
   // infeasibility measure: constraint violation
   iterate.evaluate_constraints(this->original_model);
   iterate.progress.infeasibility = this->original_model.compute_constraint_violation(iterate.evaluations.constraints, this->progress_norm);

   // optimality measure: scaled objective
   iterate.evaluate_objective(this->original_model);
   const double objective = iterate.evaluations.objective;
   iterate.progress.optimality = [=](double objective_multiplier) {
      return objective_multiplier*objective;
   };

   // auxiliary measure
   this->subproblem->set_auxiliary_measure(this->optimality_problem, iterate);
}

void FeasibilityRestoration::set_progress_measures_for_feasibility_problem(Iterate& iterate) {
   // infeasibility measure: 0
   iterate.progress.infeasibility = 0.;

   // optimality measure: constraint violation
   iterate.evaluate_constraints(this->original_model);
   const double constraint_violation = this->l1_constraint_violation_coefficient *
                                       this->original_model.compute_constraint_violation(iterate.evaluations.constraints, this->progress_norm);
   iterate.progress.optimality = [=](double /*objective_multiplier*/) {
      return constraint_violation;
   };

   // auxiliary measure
   this->subproblem->set_auxiliary_measure(this->feasibility_problem, iterate);
}

ProgressMeasures FeasibilityRestoration::compute_predicted_reduction_models_for_optimality_problem(const Iterate& current_iterate,
      const Direction& direction, double step_length) {
   // predicted infeasibility reduction: "‖c(x)‖₁ - ‖c(x) + ∇c(x)^T (αd)‖₁"
   const double current_constraint_violation = this->original_model.compute_constraint_violation(current_iterate.evaluations.constraints,
         this->progress_norm);
   const double linearized_constraint_violation = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model,
         current_iterate, direction, step_length);
   const double predicted_infeasibility_reduction = current_constraint_violation - linearized_constraint_violation;

   // predicted optimality reduction: "-∇f(x)^T (αd)"
   const double directional_derivative = dot(direction.primals, current_iterate.evaluations.objective_gradient);
   const auto predicted_optimality_reduction = [=](double objective_multiplier) {
      return step_length * (-objective_multiplier*directional_derivative);
   };

   // predicted auxiliary reduction
   const double predicted_auxiliary_reduction = this->subproblem->generate_predicted_auxiliary_reduction_model(this->optimality_problem,
         current_iterate, direction, step_length);

   return {predicted_infeasibility_reduction, predicted_optimality_reduction, predicted_auxiliary_reduction};
}

ProgressMeasures FeasibilityRestoration::compute_predicted_reduction_models_for_feasibility_problem(const Iterate& current_iterate,
      const Direction& direction, double step_length) {
   // predicted infeasibility reduction: 0
   const double predicted_infeasibility_reduction = 0.;

   // predicted optimality reduction: "‖c(x)‖₁ - ‖c(x) + ∇c(x)^T (αd)‖₁"
   const double current_constraint_violation = this->original_model.compute_constraint_violation(current_iterate.evaluations.constraints,
         this->progress_norm);
   const double linearized_constraint_violation = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model,
         current_iterate, direction, step_length);
   const auto predicted_optimality_reduction = [=](double /*objective_multiplier*/) {
      return this->l1_constraint_violation_coefficient * (current_constraint_violation - linearized_constraint_violation);
   };

   // predicted auxiliary reduction
   const double predicted_auxiliary_reduction = this->subproblem->generate_predicted_auxiliary_reduction_model(this->feasibility_problem,
         current_iterate, direction, step_length);

   return {predicted_infeasibility_reduction, predicted_optimality_reduction, predicted_auxiliary_reduction};
}

size_t FeasibilityRestoration::get_hessian_evaluation_count() const {
   return this->subproblem->get_hessian_evaluation_count();
}

size_t FeasibilityRestoration::get_number_subproblems_solved() const {
   return this->subproblem->number_subproblems_solved;
}
