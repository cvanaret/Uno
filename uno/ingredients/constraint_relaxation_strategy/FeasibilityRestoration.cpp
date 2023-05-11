// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

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
      tolerance(options.get_double("tolerance")) {
   statistics.add_column("phase", Statistics::int_width, options.get_int("statistics_restoration_phase_column_order"));
}

void FeasibilityRestoration::initialize(Iterate& initial_iterate) {
   this->subproblem->generate_initial_iterate(this->optimality_problem, initial_iterate);

   // compute the progress measures and residuals of the initial point
   this->set_infeasibility_measure(initial_iterate);
   this->set_optimality_measure(initial_iterate);
   this->subproblem->set_auxiliary_measure(this->optimality_problem, initial_iterate);
   ConstraintRelaxationStrategy::compute_primal_dual_residuals(this->optimality_problem, initial_iterate, this->residual_norm);

   // initialize the globalization strategies
   this->restoration_phase_strategy->initialize(initial_iterate);
   this->optimality_phase_strategy->initialize(initial_iterate);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
      WarmstartInformation& warmstart_information) {
   if (this->current_phase == Phase::OPTIMALITY) {
      return this->solve_optimality_problem(statistics, current_iterate, warmstart_information);
   }
   else {
      return this->solve_feasibility_problem(statistics, current_iterate, warmstart_information);
   }
}

Direction FeasibilityRestoration::solve_optimality_problem(Statistics& statistics, Iterate& current_iterate,
      WarmstartInformation& warmstart_information) {
   if (this->switched_to_optimality_phase) {
      this->switched_to_optimality_phase = false;
      warmstart_information.objective_changed = true;
      warmstart_information.constraints_changed = true;
      warmstart_information.constraint_bounds_changed = true;
      warmstart_information.variable_bounds_changed = true;
      warmstart_information.problem_changed = true;
   }

   DEBUG << "Solving the optimality subproblem\n";
   Direction direction = this->subproblem->solve(statistics, this->optimality_problem, current_iterate, warmstart_information);
   direction.objective_multiplier = 1.;
   direction.norm = norm_inf(direction.primals, Range(this->optimality_problem.number_variables));
   DEBUG2 << direction << '\n';

   // infeasible subproblem: try to minimize constraint violation by solving the feasibility subproblem
   if (direction.status == SubproblemStatus::INFEASIBLE) {
      direction = this->solve_feasibility_problem(statistics, current_iterate, direction.primals, warmstart_information);
   }
   return direction;
}

// form and solve the feasibility problem
Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
      WarmstartInformation& warmstart_information) {
   if (this->current_phase == Phase::OPTIMALITY) {
      this->switch_to_feasibility_restoration(current_iterate, warmstart_information);
   }

   DEBUG << "Solving the feasibility subproblem\n";
   Direction direction = this->subproblem->solve(statistics, this->feasibility_problem, current_iterate, warmstart_information);
   direction.objective_multiplier = 0.;
   direction.norm = norm_inf(direction.primals, Range(this->optimality_problem.number_variables));
   DEBUG2 << direction << '\n';
   assert(direction.status == SubproblemStatus::OPTIMAL && "The feasibility subproblem was not solved to optimality");
   return direction;
}

// form and solve the feasibility problem (with an initial point)
Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate, const std::vector<double>& initial_point,
      WarmstartInformation& warmstart_information) {
   this->subproblem->set_initial_point(initial_point);
   return this->solve_feasibility_problem(statistics, current_iterate, warmstart_information);
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
   if (this->current_phase == Phase::FEASIBILITY_RESTORATION &&
         ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model, current_iterate, direction,
               step_length) <= this->tolerance) {
      // if the trial infeasibility improves upon the best known infeasibility of the globalization strategy
      trial_iterate.evaluate_constraints(this->original_model);
      const double trial_infeasibility = this->original_model.compute_constraint_violation(trial_iterate.evaluations.constraints,
            this->progress_norm);
      if (this->optimality_phase_strategy->is_infeasibility_acceptable(trial_infeasibility)) {
         this->switch_to_optimality(current_iterate, trial_iterate);
      }
   }

   // evaluate the progress measures of the trial iterate
   this->set_infeasibility_measure(trial_iterate);
   this->set_optimality_measure(trial_iterate);
   this->subproblem->set_auxiliary_measure(this->current_problem(), trial_iterate);
}

void FeasibilityRestoration::switch_to_feasibility_restoration(Iterate& current_iterate, WarmstartInformation& warmstart_information) {
   DEBUG << "Switching from optimality to restoration phase\n";
   this->current_phase = Phase::FEASIBILITY_RESTORATION;
   this->optimality_phase_strategy->register_current_progress(current_iterate.progress);
   this->subproblem->initialize_feasibility_problem();
   this->subproblem->set_elastic_variable_values(this->feasibility_problem, current_iterate);
   DEBUG2 << "Current iterate:\n" << current_iterate << '\n';

   // refresh the progress measures of the current iterate
   this->set_infeasibility_measure(current_iterate);
   this->set_optimality_measure(current_iterate);
   this->subproblem->set_auxiliary_measure(this->feasibility_problem, current_iterate);

   current_iterate.multipliers.objective = 0.;
   this->restoration_phase_strategy->reset();
   this->restoration_phase_strategy->register_current_progress(current_iterate.progress);

   warmstart_information.objective_changed = true;
   warmstart_information.constraints_changed = true;
   warmstart_information.constraint_bounds_changed = true;
   warmstart_information.variable_bounds_changed = true;
   warmstart_information.problem_changed = true;
}

void FeasibilityRestoration::switch_to_optimality(Iterate& current_iterate, Iterate& trial_iterate) {
   DEBUG << "Switching from restoration to optimality phase\n";
   this->current_phase = Phase::OPTIMALITY;
   current_iterate.set_number_variables(this->optimality_problem.number_variables);
   trial_iterate.set_number_variables(this->optimality_problem.number_variables);
   this->subproblem->exit_feasibility_problem(this->optimality_problem, trial_iterate);
   this->switched_to_optimality_phase = true;

   // refresh the progress measures of current iterate
   this->set_optimality_measure(current_iterate);
   this->set_infeasibility_measure(current_iterate);
   current_iterate.multipliers.objective = 1.;
   trial_iterate.multipliers.objective = 1.;
}

bool FeasibilityRestoration::is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      double step_length) {
   // post-process the trial iterate
   this->subproblem->postprocess_iterate(this->current_problem(), trial_iterate);

   bool accept = false;
   if (direction.norm == 0.) {
      DEBUG << "Zero step acceptable\n";
      trial_iterate.evaluate_objective(this->original_model);
      accept = true;
   }
   else {
      // evaluate the predicted reduction
      ProgressMeasures predicted_reduction = {
            this->generate_predicted_infeasibility_reduction_model(current_iterate, direction, step_length),
            this->generate_predicted_optimality_reduction_model(current_iterate, direction, step_length),
            this->subproblem->generate_predicted_auxiliary_reduction_model(this->current_problem(), current_iterate, direction, step_length)
      };
      // invoke the globalization strategy for acceptance
      GlobalizationStrategy& current_phase_strategy = this->current_globalization_strategy();
      accept = current_phase_strategy.is_iterate_acceptable(statistics, trial_iterate, current_iterate.progress, trial_iterate.progress,
            predicted_reduction, this->current_problem().get_objective_multiplier());
   }

   if (accept) {
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
   return accept;
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

void FeasibilityRestoration::set_infeasibility_measure(Iterate& iterate) {
   if (this->current_phase == Phase::OPTIMALITY) {
      // constraint violation
      iterate.evaluate_constraints(this->original_model);
      iterate.progress.infeasibility = this->original_model.compute_constraint_violation(iterate.evaluations.constraints, this->progress_norm);
   }
   else {
      // 0
      iterate.progress.infeasibility = 0.;
   }
}

double FeasibilityRestoration::generate_predicted_infeasibility_reduction_model(const Iterate& current_iterate, const Direction& direction,
      double step_length) const {
   if (this->current_phase == Phase::OPTIMALITY) {
      const double current_constraint_violation = this->original_model.compute_constraint_violation(current_iterate.evaluations.constraints,
            this->progress_norm);
      const double linearized_constraint_violation = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model,
            current_iterate, direction, step_length);
      return current_constraint_violation - linearized_constraint_violation;
      //}, "‖c(x)‖₁ - ‖c(x) + ∇c(x)^T (αd)‖₁"};
   }
   else {
      return 0.;
      //}, "0"};
   }
}

void FeasibilityRestoration::set_optimality_measure(Iterate& iterate) {
   if (this->current_phase == Phase::OPTIMALITY) {
      // scaled objective
      iterate.evaluate_objective(this->original_model);
      const double objective = iterate.evaluations.objective;
      iterate.progress.optimality = [=](double objective_multiplier) {
         return objective_multiplier*objective;
      };
   }
   else {
      // constraint violation
      iterate.evaluate_constraints(this->original_model);
      const double constraint_violation = this->l1_constraint_violation_coefficient *
            this->original_model.compute_constraint_violation(iterate.evaluations.constraints, this->progress_norm);
      iterate.progress.optimality = [=](double /*objective_multiplier*/) {
         return constraint_violation;
      };
   }
}

std::function<double (double)> FeasibilityRestoration::generate_predicted_optimality_reduction_model(const Iterate& current_iterate,
      const Direction& direction, double step_length) const {
   if (this->current_phase == Phase::OPTIMALITY) {
      // precompute expensive quantities
      const double directional_derivative = dot(direction.primals, current_iterate.evaluations.objective_gradient);
      return [=](double objective_multiplier) {
         return step_length * (-objective_multiplier*directional_derivative);
      };
      //}, "-∇f(x)^T (αd)"};
   }
   else {
      const double current_constraint_violation = this->original_model.compute_constraint_violation(current_iterate.evaluations.constraints,
            this->progress_norm);
      const double linearized_constraint_violation = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model,
            current_iterate, direction, step_length);
      return [=](double /*objective_multiplier*/) {
         return this->l1_constraint_violation_coefficient * (current_constraint_violation - linearized_constraint_violation);
      };
      //}, "‖c(x)‖₁ - ‖c(x) + ∇c(x)^T (αd)‖₁"};
   }
}

size_t FeasibilityRestoration::get_hessian_evaluation_count() const {
   return this->subproblem->get_hessian_evaluation_count();
}

size_t FeasibilityRestoration::get_number_subproblems_solved() const {
   return this->subproblem->number_subproblems_solved;
}
