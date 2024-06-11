// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Options.hpp"

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
   this->compute_primal_dual_residuals(initial_iterate);
   this->set_statistics(statistics, initial_iterate);
   if (not is_linearly_feasible) {
      this->switch_to_feasibility_problem(statistics, initial_iterate);
      statistics.set("phase", "OPT");
      statistics.set("status", "linearly infeas.");
   }
   this->globalization_strategy->initialize(statistics, initial_iterate, options);
}

void FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
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
   direction.reset();
   // if we are in the optimality phase, solve the optimality problem
   if (this->current_phase == Phase::OPTIMALITY) {
      statistics.set("phase", "OPT");
      try {
         DEBUG << "Solving the optimality subproblem\n";
         this->solve_subproblem(statistics, this->optimality_problem, current_iterate, current_iterate.multipliers, direction, warmstart_information);
         if (direction.status == SubproblemStatus::INFEASIBLE) {
            // switch to the feasibility problem, starting from the current direction
            statistics.set("status", "infeas. subproblem");
            DEBUG << "/!\\ The subproblem is infeasible\n";
            this->switch_to_feasibility_problem(statistics, current_iterate);
            warmstart_information.set_cold_start();
            this->subproblem->set_initial_point(direction.primals);
         }
         else {
            return;
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
   this->solve_subproblem(statistics, this->feasibility_problem, current_iterate, current_iterate.feasibility_multipliers, direction,
         warmstart_information);
   std::swap(direction.multipliers, direction.feasibility_multipliers);
}

// an initial point is provided
void FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
      const Vector<double>& initial_point, WarmstartInformation& warmstart_information) {
   this->subproblem->set_initial_point(initial_point);
   this->compute_feasible_direction(statistics, current_iterate, direction, warmstart_information);
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

   current_iterate.set_number_variables(this->feasibility_problem.number_variables);
   this->subproblem->set_elastic_variable_values(this->feasibility_problem, current_iterate);
   DEBUG2 << "Current iterate:\n" << current_iterate << '\n';

   if (Logger::level == INFO) statistics.print_current_line();
   statistics.start_new_line();
}

void FeasibilityRestoration::solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
      const Multipliers& current_multipliers, Direction& direction, WarmstartInformation& warmstart_information) {
   // upon switching to the optimality phase, set a cold start in the subproblem solver
   if (this->switching_to_optimality_phase) {
      this->switching_to_optimality_phase = false;
      warmstart_information.set_cold_start();
   }

   direction.set_dimensions(problem.number_variables, problem.number_constraints);
   this->subproblem->solve(statistics, problem, current_iterate, current_multipliers, direction, warmstart_information);
   direction.norm = norm_inf(view(direction.primals, 0, this->model.number_variables));
   DEBUG3 << direction << '\n';
}

void FeasibilityRestoration::compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      double step_length) {
   if (this->subproblem->subproblem_definition_changed) {
      this->subproblem->subproblem_definition_changed = false;
      DEBUG << "The subproblem definition changed, the globalization strategy is reset and the auxiliary measure is recomputed\n";
      this->globalization_strategy->reset();
      this->subproblem->set_auxiliary_measure(this->model, current_iterate);
   }
   this->evaluate_progress_measures(trial_iterate);

   // possibly go from restoration phase to optimality phase
   if (this->current_phase == Phase::FEASIBILITY_RESTORATION && not this->switch_to_optimality_requires_acceptance &&
         (not this->switch_to_optimality_requires_linearized_feasibility || this->model.constraint_violation(current_iterate.evaluations.constraints +
         step_length*(current_iterate.evaluations.constraint_jacobian * direction.primals), this->residual_norm) <= this->linear_feasibility_tolerance) &&
         this->globalization_strategy->is_infeasibility_sufficiently_reduced(current_iterate.progress, trial_iterate.progress)) {
      this->switch_to_optimality_phase(current_iterate, trial_iterate);
   }
}

void FeasibilityRestoration::switch_to_optimality_phase(Iterate& current_iterate, Iterate& trial_iterate) {
   DEBUG << "Switching from restoration to optimality phase\n";
   this->current_phase = Phase::OPTIMALITY;
   this->globalization_strategy->register_current_progress(current_iterate.progress);
   current_iterate.set_number_variables(this->optimality_problem.number_variables);
   trial_iterate.set_number_variables(this->optimality_problem.number_variables);
   current_iterate.objective_multiplier = trial_iterate.objective_multiplier = 1.;

   this->subproblem->exit_feasibility_problem(this->optimality_problem, trial_iterate);
   this->switching_to_optimality_phase = true;
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
      DEBUG << "Zero step acceptable\n";
      trial_iterate.evaluate_objective(this->model);
      accept_iterate = true;
      statistics.set("status", "accepted (0 primal step)");
   }
   else {
      // invoke the globalization strategy for acceptance
      const ProgressMeasures predicted_reduction = this->compute_predicted_reduction_models(current_iterate, direction, step_length);
      accept_iterate = this->globalization_strategy->is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
            predicted_reduction, this->current_problem().get_objective_multiplier());
   }
   // possibly switch back to optimality phase if the trial iterate in feasibility restoration was accepted
   if (accept_iterate) {
      if (this->current_phase == Phase::FEASIBILITY_RESTORATION && this->switch_to_optimality_requires_acceptance &&
            this->globalization_strategy->is_infeasibility_sufficiently_reduced(current_iterate.progress, trial_iterate.progress)) {
         this->switch_to_optimality_phase(current_iterate, trial_iterate);
      }
      this->compute_primal_dual_residuals(trial_iterate);
      this->set_dual_residuals_statistics(statistics, trial_iterate);
   }
   ConstraintRelaxationStrategy::set_progress_statistics(statistics, trial_iterate);
   return accept_iterate;
}

void FeasibilityRestoration::compute_primal_dual_residuals(Iterate& iterate) {
   ConstraintRelaxationStrategy::compute_primal_dual_residuals(this->optimality_problem, this->feasibility_problem, iterate);
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
   this->subproblem->set_auxiliary_measure(this->model, iterate);
}

ProgressMeasures FeasibilityRestoration::compute_predicted_reduction_models(Iterate& current_iterate, const Direction& direction, double step_length) {
   return {
      this->compute_predicted_infeasibility_reduction_model(current_iterate, direction.primals, step_length),
      this->compute_predicted_objective_reduction_model(current_iterate, direction.primals, step_length, this->subproblem->get_lagrangian_hessian()),
      this->subproblem->compute_predicted_auxiliary_reduction_model(this->model, current_iterate, direction.primals, step_length)
   };
}

void FeasibilityRestoration::set_trust_region_radius(double trust_region_radius) {
   this->subproblem->set_trust_region_radius(trust_region_radius);
}

size_t FeasibilityRestoration::maximum_number_variables() const {
   return std::max(this->optimality_problem.number_variables, this->feasibility_problem.number_variables);
}

size_t FeasibilityRestoration::maximum_number_constraints() const {
   return std::max(this->optimality_problem.number_constraints, this->feasibility_problem.number_constraints);
}

void FeasibilityRestoration::set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const {
   if (this->current_phase == Phase::OPTIMALITY) {
      statistics.set("complementarity", iterate.residuals.complementarity);
      statistics.set("stationarity", iterate.residuals.KKT_stationarity);
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
