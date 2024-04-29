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
      subproblem(SubproblemFactory::create(
         // allocate the largest size necessary to solve the optimality subproblem or the feasibility subproblem
         std::max(this->optimality_problem.number_variables, this->feasibility_problem.number_variables),
         std::max(this->optimality_problem.number_constraints, this->feasibility_problem.number_constraints),
         std::max(this->optimality_problem.number_objective_gradient_nonzeros(), this->feasibility_problem.number_objective_gradient_nonzeros()),
         std::max(this->optimality_problem.number_jacobian_nonzeros(), this->feasibility_problem.number_jacobian_nonzeros()),
         std::max(this->optimality_problem.number_hessian_nonzeros(), this->feasibility_problem.number_hessian_nonzeros()),
         options)),
      globalization_strategy(GlobalizationStrategyFactory::create(options.get_string("globalization_strategy"), options)),
      linear_feasibility_tolerance(options.get_double("tolerance")),
      switch_to_optimality_requires_acceptance(options.get_bool("switch_to_optimality_requires_acceptance")),
      switch_to_optimality_requires_linearized_feasibility(options.get_bool("switch_to_optimality_requires_linearized_feasibility")) {
}

void FeasibilityRestoration::initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) {
   // statistics
   this->subproblem->initialize_statistics(statistics, options);
   statistics.add_column("phase", Statistics::int_width, options.get_int("statistics_restoration_phase_column_order"));
   statistics.set("phase", "OPT");

   // initial iterate
   const bool is_linearly_feasible = this->subproblem->generate_initial_iterate(this->optimality_problem, initial_iterate);
   this->evaluate_progress_measures(initial_iterate);
   this->compute_primal_dual_residuals(this->feasibility_problem, initial_iterate);
   this->set_statistics(statistics, initial_iterate);
   if (not is_linearly_feasible) {
      this->switch_to_feasibility_problem(statistics, initial_iterate);
      statistics.set("phase", "OPT");
      statistics.set("status", "linearly infeas.");
   }
   this->globalization_strategy->initialize(statistics, initial_iterate, options);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
      WarmstartInformation& warmstart_information) {
   /*
   if (1e6 < norm_inf(current_iterate.multipliers.constraints)) {
      // large duals are an indication of CQ failure
      statistics.set("status", "large duals");
      DEBUG << "/!\\ The duals are large\n";
      this->switch_to_feasibility_problem(statistics, current_iterate);
      warmstart_information.set_cold_start();
   }
   */
   // if we are in the optimality phase, solve the optimality problem
   if (this->current_phase == Phase::OPTIMALITY) {
      statistics.set("phase", "OPT");
      try {
         DEBUG << "Solving the optimality subproblem\n";
         Direction direction = this->solve_subproblem(statistics, this->optimality_problem, current_iterate, warmstart_information);
         if (direction.status == SubproblemStatus::INFEASIBLE) {
            // switch to the feasibility problem, starting from the current direction
            statistics.set("status", "infeas. subproblem");
            DEBUG << "/!\\ The subproblem is infeasible\n";
            this->switch_to_feasibility_problem(statistics, current_iterate);
            warmstart_information.set_cold_start();
            this->subproblem->set_initial_point(direction.primals);
         }
         else {
            return direction;
         }
      }
      catch (const UnstableRegularization&) {
         this->switch_to_feasibility_problem(statistics, current_iterate);
         warmstart_information.set_cold_start();
      }
   }

   // solve the feasibility problem (minimize the constraint violation)
   DEBUG << "Solving the feasibility subproblem\n";
   statistics.set("phase", "FEAS");
   // note: failure of regularization should not happen here, since the feasibility Jacobian has full rank
   return this->solve_subproblem(statistics, this->feasibility_problem, current_iterate, warmstart_information);
}

// an initial point is provided
Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
      const std::vector<double>& initial_point, WarmstartInformation& warmstart_information) {
   this->subproblem->set_initial_point(initial_point);
   return this->compute_feasible_direction(statistics, current_iterate, warmstart_information);
}

bool FeasibilityRestoration::solving_feasibility_problem() const {
   return (this->current_phase == Phase::FEASIBILITY_RESTORATION);
}

// precondition: this->current_phase == Phase::OPTIMALITY
void FeasibilityRestoration::switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate) {
   DEBUG << "Switching from optimality to restoration phase\n";
   this->current_phase = Phase::FEASIBILITY_RESTORATION;
   this->globalization_strategy->register_current_progress(current_iterate.progress);
   this->subproblem->initialize_feasibility_problem(this->feasibility_problem, current_iterate);

   // reset multipliers (this means no curvature in the first feasibility restoration iteration)
   current_iterate.multipliers.reset();
   // TODO: allocate the iterates with the maximum number of dimensions and never physically resize
   current_iterate.set_number_variables(this->feasibility_problem.number_variables);
   this->subproblem->set_elastic_variable_values(this->feasibility_problem, current_iterate);
   DEBUG2 << "Current iterate:\n" << current_iterate << '\n';

   if (Logger::level == INFO) statistics.print_current_line();
   statistics.start_new_line();
}

Direction FeasibilityRestoration::solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
      WarmstartInformation& warmstart_information) {
   // upon switching to the optimality phase, set a cold start in the subproblem solver
   if (this->switching_to_optimality_phase) {
      this->switching_to_optimality_phase = false;
      warmstart_information.set_cold_start();
   }

   Direction direction = this->subproblem->solve(statistics, problem, current_iterate, warmstart_information);
   direction.norm = norm_inf(view(direction.primals, this->model.number_variables));
   direction.multipliers.objective = problem.get_objective_multiplier();
   DEBUG3 << direction << '\n';
   return direction;
}

void FeasibilityRestoration::compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      double step_length) {
   if (this->subproblem->subproblem_definition_changed) {
      this->subproblem->subproblem_definition_changed = false;
      DEBUG << "The subproblem definition changed, the globalization strategy is reset and the auxiliary measure is recomputed\n";
      this->globalization_strategy->reset();
      this->subproblem->set_auxiliary_measure(this->current_problem(), current_iterate);
   }
   this->evaluate_progress_measures(trial_iterate);

   // possibly go from restoration phase to optimality phase
   if (this->current_phase == Phase::FEASIBILITY_RESTORATION && not this->switch_to_optimality_requires_acceptance &&
         (not this->switch_to_optimality_requires_linearized_feasibility || this->model.linearized_constraint_violation(direction.primals,
         current_iterate.evaluations.constraints, current_iterate.evaluations.constraint_jacobian, step_length, this->residual_norm) <=
         this->linear_feasibility_tolerance) && this->globalization_strategy->is_feasibility_iterate_acceptable(current_iterate.progress,
               trial_iterate.progress)) {
      this->switch_to_optimality_phase(current_iterate, trial_iterate);
   }
}

void FeasibilityRestoration::switch_to_optimality_phase(Iterate& current_iterate, Iterate& trial_iterate) {
   DEBUG << "Switching from restoration to optimality phase\n";
   this->current_phase = Phase::OPTIMALITY;
   current_iterate.set_number_variables(this->optimality_problem.number_variables);
   trial_iterate.set_number_variables(this->optimality_problem.number_variables);
   this->subproblem->exit_feasibility_problem(this->optimality_problem, trial_iterate);
   this->switching_to_optimality_phase = true;
   this->globalization_strategy->register_current_progress(current_iterate.progress);

   current_iterate.multipliers.objective = trial_iterate.multipliers.objective = FeasibilityRestoration::objective_multiplier;
}

bool FeasibilityRestoration::is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      double step_length) {
   this->subproblem->postprocess_iterate(this->current_problem(), trial_iterate);
   this->compute_progress_measures(current_iterate, trial_iterate, direction, step_length);

   bool accept_iterate = false;
   /*
   // detect trivial multipliers
   if (this->current_phase == Phase::FEASIBILITY_RESTORATION && not trial_iterate.multipliers.not_all_zero(this->model.number_variables, 1e-6)) {
      DEBUG << "Trivial duals in restoration phase, trial iterate rejected\n\n";
      statistics.set("status", "rejected (trivial duals)");
   }
   else
   */
   if (direction.norm == 0.) {
      DEBUG << "Zero step acceptable\n\n";
      trial_iterate.evaluate_objective(this->model);
      accept_iterate = true;
      statistics.set("status", "accepted (0 step)");
   }
   else {
      // invoke the globalization strategy for acceptance
      const ProgressMeasures predicted_reduction = this->compute_predicted_reduction_models(current_iterate, direction, step_length);
      accept_iterate = this->globalization_strategy->is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
            predicted_reduction, this->current_problem().get_objective_multiplier());
   }
   // possibly switch back to optimality phase if the trial iterate in feasibility restoration was accepted
   if (accept_iterate && this->current_phase == Phase::FEASIBILITY_RESTORATION && this->switch_to_optimality_requires_acceptance &&
         this->globalization_strategy->is_feasibility_iterate_acceptable(current_iterate.progress, trial_iterate.progress)) {
      this->switch_to_optimality_phase(current_iterate, trial_iterate);
   }
   this->compute_primal_dual_residuals(this->feasibility_problem, trial_iterate);
   this->set_statistics(statistics, trial_iterate);
   return accept_iterate;
}

const OptimizationProblem& FeasibilityRestoration::current_problem() const {
   if (this->current_phase == Phase::OPTIMALITY) {
      return this->optimality_problem;
   }
   else {
      return this->feasibility_problem;
   }
}

void FeasibilityRestoration::evaluate_progress_measures(Iterate& iterate) const {
   this->set_infeasibility_measure(iterate);
   this->set_objective_measure(iterate);
   this->subproblem->set_auxiliary_measure(this->optimality_problem, iterate);
}

ProgressMeasures FeasibilityRestoration::compute_predicted_reduction_models(Iterate& current_iterate, const Direction& direction, double step_length) {
   return {
      this->compute_predicted_infeasibility_reduction_model(current_iterate, direction, step_length),
      this->compute_predicted_objective_reduction_model(current_iterate, direction, step_length, this->subproblem->get_lagrangian_hessian()),
      this->subproblem->compute_predicted_auxiliary_reduction_model(this->optimality_problem, current_iterate, direction, step_length)
   };
}

void FeasibilityRestoration::set_trust_region_radius(double trust_region_radius) {
   this->subproblem->set_trust_region_radius(trust_region_radius);
}

double FeasibilityRestoration::complementarity_error(const std::vector<double>& primals, const std::vector<double>& constraints,
      const Multipliers& multipliers) const {
   // bound constraints
   const VectorExpression<double, Range<FORWARD>> variable_complementarity(Range(this->model.number_variables), [&](size_t variable_index) {
      if (0. < multipliers.lower_bounds[variable_index]) {
         return multipliers.lower_bounds[variable_index] * (primals[variable_index] - this->model.variable_lower_bound(variable_index));
      }
      if (multipliers.upper_bounds[variable_index] < 0.) {
         return multipliers.upper_bounds[variable_index] * (primals[variable_index] - this->model.variable_upper_bound(variable_index));
      }
      return 0.;
   });

   // inequality constraints
   const VectorExpression<double, const Collection<size_t>&> constraint_complementarity(this->model.get_inequality_constraints(), [&](size_t
   constraint_index) {
      if (0. < multipliers.constraints[constraint_index]) { // lower bound
         return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->model.constraint_lower_bound(constraint_index));
      }
      else if (multipliers.constraints[constraint_index] < 0.) { // upper bound
         return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->model.constraint_upper_bound(constraint_index));
      }
      return 0.;
   });
   return norm(this->residual_norm, variable_complementarity, constraint_complementarity);
}

void FeasibilityRestoration::set_statistics(Statistics& statistics, const Iterate& iterate) const {
   statistics.set("objective", iterate.evaluations.objective);
   if (this->model.is_constrained()) {
      statistics.set("primal infeas.", iterate.residuals.infeasibility);
   }
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
