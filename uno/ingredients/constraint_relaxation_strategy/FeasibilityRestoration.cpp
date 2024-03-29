// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "linear_algebra/view.hpp"

FeasibilityRestoration::FeasibilityRestoration(const Model& model, const Options& options) :
      ConstraintRelaxationStrategy(model, options),
      // create the (optimality phase) optimality problem (= original model)
      optimality_problem(model),
      // create the (restoration phase) feasibility problem (objective multiplier = 0)
      feasibility_problem(model, 0., options.get_double("l1_constraint_violation_coefficient")),
      subproblem(SubproblemFactory::create(this->feasibility_problem.number_variables, this->feasibility_problem.number_constraints,
            this->feasibility_problem.number_objective_gradient_nonzeros(), this->feasibility_problem.number_jacobian_nonzeros(),
            this->feasibility_problem.number_hessian_nonzeros(), options)),
      // create the globalization strategies (one for each phase)
      restoration_phase_strategy(GlobalizationStrategyFactory::create(options.get_string("globalization_strategy"), false, options)),
      optimality_phase_strategy(GlobalizationStrategyFactory::create(options.get_string("globalization_strategy"), true, options)),
      tolerance(options.get_double("tolerance")),
      switch_to_optimality_requires_acceptance(options.get_bool("switch_to_optimality_requires_acceptance")),
      switch_to_optimality_requires_linearized_feasibility(options.get_bool("switch_to_optimality_requires_linearized_feasibility")) {
}

void FeasibilityRestoration::initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) {
   // compute the progress measures and residuals of the initial point
   this->subproblem->generate_initial_iterate(this->optimality_problem, initial_iterate);
   this->evaluate_progress_measures(this->optimality_problem, initial_iterate);
   this->compute_primal_dual_residuals(this->original_model, this->feasibility_problem, initial_iterate);

   // initialize the globalization strategies
   this->restoration_phase_strategy->initialize(statistics, initial_iterate, options);
   this->optimality_phase_strategy->initialize(statistics, initial_iterate, options);

   this->subproblem->initialize_statistics(statistics, options);
   this->set_statistics(statistics, initial_iterate);
   statistics.add_column("phase", Statistics::int_width, options.get_int("statistics_restoration_phase_column_order"));
   statistics.set("phase", static_cast<int>(this->current_phase));
   if (this->original_model.is_constrained()) {
      statistics.set("primal infeas.", initial_iterate.progress.infeasibility);
   }
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
      WarmstartInformation& warmstart_information) {
   // if we are in the optimality phase, solve the optimality problem
   if (this->current_phase == Phase::OPTIMALITY) {
      statistics.set("phase", static_cast<int>(this->current_phase));
      try {
         DEBUG << "Solving the optimality subproblem\n";
         Direction direction = this->solve_subproblem(statistics, this->optimality_problem, current_iterate, warmstart_information);
         if (direction.status != SubproblemStatus::INFEASIBLE) {
            return direction;
         }
         else {
            // infeasible subproblem: switch to the feasibility problem, starting from the current direction
            statistics.set("status", "infeas. subproblem");
            DEBUG << "/!\\ The subproblem is infeasible\n";
            this->switch_to_feasibility_problem(statistics, current_iterate, warmstart_information);
            this->subproblem->set_initial_point(direction.primals);
         }
      }
      catch (const UnstableRegularization&) {
         this->switch_to_feasibility_problem(statistics, current_iterate, warmstart_information);
      }
   }

   // solve the feasibility problem (minimize the constraint violation)
   DEBUG << "Solving the feasibility subproblem\n";
   statistics.set("phase", static_cast<int>(this->current_phase));
   // note: failure of regularization should not happen here, since the feasibility Jacobian has full rank
   return this->solve_subproblem(statistics, this->feasibility_problem, current_iterate, warmstart_information);
}

// an initial point is provided
Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
      const std::vector<double>& initial_point, WarmstartInformation& warmstart_information) {
   this->subproblem->set_initial_point(initial_point);
   return this->compute_feasible_direction(statistics, current_iterate, warmstart_information);
}

bool FeasibilityRestoration::solving_feasibility_problem() {
   return (this->current_phase == Phase::FEASIBILITY_RESTORATION);
}

void FeasibilityRestoration::switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate, WarmstartInformation& warmstart_information) {
   if (this->current_phase == Phase::FEASIBILITY_RESTORATION) {
      throw std::runtime_error("FeasibilityRestoration::switch_to_feasibility_problem: already in feasibility restoration.\n");
   }
   
   DEBUG << "Switching from optimality to restoration phase\n";
   this->current_phase = Phase::FEASIBILITY_RESTORATION;
   this->optimality_phase_strategy->register_current_progress(current_iterate.progress);
   this->subproblem->initialize_feasibility_problem();
   this->subproblem->set_elastic_variable_values(this->feasibility_problem, current_iterate);
   DEBUG2 << "Current iterate:\n" << current_iterate << '\n';

   // compute the progress measures of the current iterate for the feasibility problem
   this->evaluate_progress_measures(this->feasibility_problem, current_iterate);
   current_iterate.multipliers.objective = 0.;

   this->restoration_phase_strategy->reset();
   this->restoration_phase_strategy->register_current_progress(current_iterate.progress);
   warmstart_information.set_cold_start();
   
   if (Logger::level == INFO) statistics.print_current_line();
   statistics.start_new_line();
}

Direction FeasibilityRestoration::solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
      WarmstartInformation& warmstart_information) {
   if (this->switched_to_optimality_phase) {
      this->switched_to_optimality_phase = false;
      warmstart_information.set_cold_start();
   }

   Direction direction = this->subproblem->solve(statistics, problem, current_iterate, warmstart_information);
   direction.norm = norm_inf(view(direction.primals, this->original_model.number_variables));
   direction.multipliers.objective = problem.get_objective_multiplier();
   DEBUG3 << direction << '\n';
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
   if (this->current_phase == Phase::FEASIBILITY_RESTORATION && not this->switch_to_optimality_requires_acceptance &&
         (not this->switch_to_optimality_requires_linearized_feasibility || this->original_model.linearized_constraint_violation(direction.primals,
               current_iterate.evaluations.constraints, current_iterate.evaluations.constraint_jacobian, step_length, this->residual_norm) <= this->tolerance)) {
      // if the trial infeasibility improves upon the best known infeasibility of the globalization strategy
      trial_iterate.evaluate_constraints(this->original_model);
      const double trial_infeasibility = this->original_model.constraint_violation(trial_iterate.evaluations.constraints, this->progress_norm);
      if (this->optimality_phase_strategy->is_infeasibility_acceptable(trial_infeasibility)) {
         this->switch_to_optimality_phase(current_iterate, trial_iterate);
      }
   }

   // evaluate the progress measures of the trial iterate
   this->evaluate_progress_measures(this->current_problem(), trial_iterate);
}

void FeasibilityRestoration::switch_to_optimality_phase(Iterate& current_iterate, Iterate& trial_iterate) {
   DEBUG << "Switching from restoration to optimality phase\n";
   this->current_phase = Phase::OPTIMALITY;
   current_iterate.set_number_variables(this->optimality_problem.number_variables);
   trial_iterate.set_number_variables(this->optimality_problem.number_variables);
   this->subproblem->exit_feasibility_problem(this->optimality_problem, trial_iterate);
   this->switched_to_optimality_phase = true;

   // refresh the progress measures of current iterate
   this->evaluate_progress_measures(this->optimality_problem, current_iterate);
   current_iterate.multipliers.objective = 1.;
   trial_iterate.multipliers.objective = 1.;
}

bool FeasibilityRestoration::is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      double step_length) {
   this->subproblem->postprocess_iterate(this->current_problem(), trial_iterate);
   this->compute_progress_measures(current_iterate, trial_iterate, direction, step_length);

   bool accept_iterate = false;
   if (direction.norm == 0.) {
      DEBUG << "Zero step acceptable\n\n";
      trial_iterate.evaluate_objective(this->original_model);
      accept_iterate = true;
      statistics.set("status", "accepted (0 step)");
   }
   else {
      // evaluate the predicted reduction
      ProgressMeasures predicted_reduction = this->compute_predicted_reduction_models(current_iterate, direction, step_length);

      // invoke the globalization strategy for acceptance
      accept_iterate = this->current_globalization_strategy().is_iterate_acceptable(statistics, trial_iterate, current_iterate.progress,
            trial_iterate.progress, predicted_reduction, this->current_problem().get_objective_multiplier());
   }

   if (accept_iterate) {
      this->compute_primal_dual_residuals(this->original_model, this->feasibility_problem, trial_iterate);
      this->set_statistics(statistics, trial_iterate);
   }
   if (this->original_model.is_constrained()) {
      statistics.set("primal infeas.", trial_iterate.progress.infeasibility);
   }
   return accept_iterate;
}

void FeasibilityRestoration::evaluate_progress_measures(const OptimizationProblem& problem, Iterate& iterate) const {
   problem.set_infeasibility_measure(iterate, this->progress_norm);
   problem.set_optimality_measure(iterate);
   this->subproblem->set_auxiliary_measure(problem, iterate);
}

ProgressMeasures FeasibilityRestoration::compute_predicted_reduction_models(Iterate& current_iterate, const Direction& direction, double step_length) {
   const OptimizationProblem& current_problem = this->current_problem();
   return {
      current_problem.compute_predicted_infeasibility_reduction_model(current_iterate, direction, step_length, this->progress_norm),
      this->subproblem->compute_predicted_optimality_reduction_model(current_problem, current_iterate, direction, step_length),
      this->subproblem->compute_predicted_auxiliary_reduction_model(current_problem, current_iterate, direction, step_length)
   };
}

const OptimizationProblem& FeasibilityRestoration::current_problem() const {
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

double FeasibilityRestoration::compute_complementarity_error(const std::vector<double>& primals, const std::vector<double>& constraints,
      const Multipliers& multipliers) const {
   // bound constraints
   VectorExpression<double> variable_complementarity(this->original_model.number_variables, [&](size_t variable_index) {
      if (0. < multipliers.lower_bounds[variable_index]) {
         return multipliers.lower_bounds[variable_index] * (primals[variable_index] - this->original_model.variable_lower_bound(variable_index));
      }
      if (multipliers.upper_bounds[variable_index] < 0.) {
         return multipliers.upper_bounds[variable_index] * (primals[variable_index] - this->original_model.variable_upper_bound(variable_index));
      }
      return 0.;
   });

   // constraints
   VectorExpression<double> constraint_complementarity(this->original_model.inequality_constraints.size(), [&](size_t inequality_index) {
      const size_t constraint_index = this->original_model.inequality_constraints[inequality_index];
      if (0. < multipliers.constraints[constraint_index]) { // lower bound
         return multipliers.constraints[constraint_index] * (constraints[constraint_index] -
               this->original_model.constraint_lower_bound(constraint_index));
      }
      else if (multipliers.constraints[constraint_index] < 0.) { // upper bound
         return multipliers.constraints[constraint_index] * (constraints[constraint_index] -
               this->original_model.constraint_upper_bound(constraint_index));
      }
      return 0.;
   });
   return norm(this->residual_norm, variable_complementarity, constraint_complementarity);
}

void FeasibilityRestoration::set_statistics(Statistics& statistics, const Iterate& iterate) const {
   if (this->current_phase == Phase::OPTIMALITY) {
      statistics.set("complementarity", iterate.residuals.optimality_complementarity);
      statistics.set("stationarity", iterate.residuals.optimality_stationarity);
   }
   else {
      statistics.set("complementarity", iterate.residuals.feasibility_complementarity);
      statistics.set("stationarity", iterate.residuals.feasibility_stationarity);
   }
}

size_t FeasibilityRestoration::get_hessian_evaluation_count() const {
   return this->subproblem->get_hessian_evaluation_count();
}

size_t FeasibilityRestoration::get_number_subproblems_solved() const {
   return this->subproblem->number_subproblems_solved;
}
