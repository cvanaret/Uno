// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "l1Relaxation.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"
#include "linear_algebra/view.hpp"

/*
 * Infeasibility detection and SQP methods for nonlinear optimization
 * Richard H. Byrd, Frank E. Curtis and Jorge Nocedal
 * http://epubs.siam.org/doi/pdf/10.1137/080738222
 */

l1Relaxation::l1Relaxation(Statistics& statistics, const Model& model, const Options& options) :
      ConstraintRelaxationStrategy(model, options),
      // create the l1 feasibility problem (objective multiplier = 0)
      feasibility_problem(model, 0., options.get_double("l1_constraint_violation_coefficient")),
      // create the l1 relaxed problem
      l1_relaxed_problem(model, options.get_double("l1_relaxation_initial_parameter"), options.get_double("l1_constraint_violation_coefficient")),
      subproblem(SubproblemFactory::create(statistics, this->l1_relaxed_problem.number_variables, this->l1_relaxed_problem.number_constraints,
            this->l1_relaxed_problem.get_number_jacobian_nonzeros(), this->l1_relaxed_problem.get_number_hessian_nonzeros(), options)),
      globalization_strategy(GlobalizationStrategyFactory::create(statistics, options.get_string("globalization_strategy"), true, options)),
      penalty_parameter(options.get_double("l1_relaxation_initial_parameter")),
      tolerance(options.get_double("tolerance")),
      parameters({
         options.get_bool("l1_relaxation_fixed_parameter"),
         options.get_double("l1_relaxation_decrease_factor"),
         options.get_double("l1_relaxation_epsilon1"),
         options.get_double("l1_relaxation_epsilon2"),
         options.get_double("l1_relaxation_residual_small_threshold")
      }),
      small_duals_threshold(options.get_double("l1_small_duals_threshold")),
      trial_multipliers(this->l1_relaxed_problem.number_variables, model.number_constraints) {
   statistics.add_column("penalty param.", Statistics::double_width, options.get_int("statistics_penalty_parameter_column_order"));
}

void l1Relaxation::initialize(Iterate& initial_iterate) {
   this->subproblem->set_elastic_variable_values(this->l1_relaxed_problem, initial_iterate);
   this->subproblem->generate_initial_iterate(this->l1_relaxed_problem, initial_iterate);

   // compute the progress measures and residuals of the initial point
   this->set_progress_measures(initial_iterate);
   this->compute_primal_dual_residuals(this->original_model, this->feasibility_problem, initial_iterate);

   // initialize the globalization strategy
   this->globalization_strategy->initialize(initial_iterate);
}

Direction l1Relaxation::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, WarmstartInformation& warmstart_information) {
   if (0. < this->penalty_parameter) {
      return this->solve_sequence_of_relaxed_subproblems(statistics, current_iterate, warmstart_information);
   }
   else {
      return this->solve_subproblem(statistics, this->feasibility_problem, current_iterate, warmstart_information);
   }
}

// an initial point is provided
Direction l1Relaxation::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, const std::vector<double>& initial_point,
      WarmstartInformation& warmstart_information) {
   this->subproblem->set_initial_point(initial_point);
   return this->compute_feasible_direction(statistics, current_iterate, warmstart_information);
}

void l1Relaxation::switch_to_feasibility_problem(Iterate& /*current_iterate*/, WarmstartInformation& /*warmstart_information*/) {
   throw std::runtime_error("l1Relaxation::switch_to_feasibility_problem is not implemented\n");
}

// use Byrd's steering rules to update the penalty parameter and compute a descent direction
Direction l1Relaxation::solve_sequence_of_relaxed_subproblems(Statistics& statistics, Iterate& current_iterate,
      WarmstartInformation& warmstart_information) {
   // stage a: compute a direction for the current penalty parameter
   Direction direction = this->solve_l1_relaxed_problem(statistics, current_iterate, this->penalty_parameter, warmstart_information);
   // from now on, only the penalty parameter, therefore the objective, changes
   warmstart_information.only_objective_changed();

   // penalty update: if penalty parameter is already 0 or fixed by the user, no need to decrease it
   if (0. < this->penalty_parameter && not this->parameters.fixed_parameter) {
      double linearized_residual = this->original_model.compute_linearized_constraint_violation(direction.primals,
            current_iterate.evaluations.constraints, current_iterate.evaluations.constraint_jacobian, direction.primal_dual_step_length, Norm::L1);
      DEBUG << "Linearized infeasibility mk(dk): " << linearized_residual << "\n\n";

      // if the current direction is already feasible, terminate
      if (this->tolerance < linearized_residual) {
         const double current_penalty_parameter = this->penalty_parameter;

         // stage c: compute the lowest possible constraint violation (penalty parameter = 0)
         DEBUG << "Compute ideal solution by solving the feasibility problem:\n";
         this->subproblem->initialize_feasibility_problem();
         Direction feasibility_direction = this->solve_subproblem(statistics, this->feasibility_problem, current_iterate, warmstart_information);
         const double residual_lowest_violation = this->original_model.compute_linearized_constraint_violation(feasibility_direction.primals,
               current_iterate.evaluations.constraints, current_iterate.evaluations.constraint_jacobian,
               feasibility_direction.primal_dual_step_length, Norm::L1);
         DEBUG << "Lowest linearized infeasibility mk(dk): " << residual_lowest_violation << '\n';
         // TODO let the subproblem exit the feasibility problem

         // stage f: update the penalty parameter based on the current dual error
         this->decrease_parameter_aggressively(current_iterate, feasibility_direction);
         if (this->penalty_parameter == 0.) {
            direction = feasibility_direction;
         }
         else {
            if (this->penalty_parameter < current_penalty_parameter) {
               direction = this->solve_l1_relaxed_problem(statistics, current_iterate, this->penalty_parameter, warmstart_information);
               linearized_residual = this->original_model.compute_linearized_constraint_violation(direction.primals,
                     current_iterate.evaluations.constraints, current_iterate.evaluations.constraint_jacobian, direction.primal_dual_step_length,
                     Norm::L1);
            }

            // stage d: further decrease penalty parameter to reach a fraction of the ideal decrease
            direction = this->enforce_linearized_residual_sufficient_decrease(statistics, current_iterate, direction, linearized_residual,
                  residual_lowest_violation, warmstart_information);
            // stage e: further decrease penalty parameter to guarantee a descent direction for the l1 merit function
            direction = this->enforce_descent_direction_for_l1_merit(statistics, current_iterate, direction, feasibility_direction,
                  warmstart_information);
         }
      }
   }
   return direction;
}

Direction l1Relaxation::solve_subproblem(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
      const WarmstartInformation& warmstart_information) {
   DEBUG << "Solving the subproblem with penalty parameter " << problem.get_objective_multiplier() << "\n\n";

   // solve the subproblem
   Direction direction = this->subproblem->solve(statistics, problem, current_iterate, warmstart_information);
   direction.norm = norm_inf(view(direction.primals, this->original_model.number_variables));
   direction.multipliers.objective = problem.get_objective_multiplier();
   DEBUG2 << direction << '\n';
   assert(direction.status == SubproblemStatus::OPTIMAL && "The subproblem was not solved to optimality");
   return direction;
}

Direction l1Relaxation::solve_l1_relaxed_problem(Statistics& statistics, Iterate& current_iterate, double current_penalty_parameter,
      const WarmstartInformation& warmstart_information) {
   this->l1_relaxed_problem.set_objective_multiplier(current_penalty_parameter);
   return this->solve_subproblem(statistics, this->l1_relaxed_problem, current_iterate, warmstart_information);
}

void l1Relaxation::decrease_parameter_aggressively(Iterate& current_iterate, const Direction& direction) {
   add_vectors(current_iterate.multipliers.constraints, direction.multipliers.constraints, direction.primal_dual_step_length,
         this->trial_multipliers.constraints);
   add_vectors(current_iterate.multipliers.lower_bounds, direction.multipliers.lower_bounds, direction.bound_dual_step_length,
         this->trial_multipliers.lower_bounds);
   add_vectors(current_iterate.multipliers.upper_bounds, direction.multipliers.upper_bounds, direction.bound_dual_step_length,
         this->trial_multipliers.upper_bounds);

   // there must be at least a nonzero dual to avoid trivial stationary points
   if (this->trial_multipliers.not_all_zero(this->original_model.number_variables, this->small_duals_threshold)) {
      // compute the ideal error (with a zero penalty parameter)
      const double infeasible_dual_error = l1Relaxation::compute_infeasible_dual_error(current_iterate);
      DEBUG << "Ideal dual error: " << infeasible_dual_error << '\n';
      const double scaled_error = infeasible_dual_error / std::max(1., current_iterate.residuals.infeasibility);
      this->penalty_parameter = std::min(this->penalty_parameter, scaled_error * scaled_error);
      DEBUG << "Further aggressively decrease the penalty parameter to " << this->penalty_parameter << '\n';
   }
   else {
      WARNING << RED << "l1Relaxation: all multipliers are almost 0. The penalty parameter won't be decreased" << RESET << '\n';
   }
}

// measure that combines KKT error and complementarity error
double l1Relaxation::compute_infeasible_dual_error(Iterate& current_iterate) {
   // stationarity error
   ConstraintRelaxationStrategy::evaluate_lagrangian_gradient(this->original_model.number_variables, current_iterate, this->trial_multipliers, 0.);
   double error = norm_1(current_iterate.lagrangian_gradient.constraints_contribution);

   // complementarity error
   error += this->feasibility_problem.compute_complementarity_error(current_iterate.primals, current_iterate.evaluations.constraints,
         this->trial_multipliers, Norm::L1);
   return error;
}

Direction l1Relaxation::enforce_linearized_residual_sufficient_decrease(Statistics& statistics, Iterate& current_iterate, Direction& direction,
      double linearized_residual, double residual_lowest_violation, WarmstartInformation& warmstart_information) {
   while (0. < this->penalty_parameter && not this->linearized_residual_sufficient_decrease(current_iterate, linearized_residual,
         residual_lowest_violation)) {
      // decrease the penalty parameter and re-solve the problem
      this->penalty_parameter /= this->parameters.decrease_factor;
      DEBUG << "Further decrease the penalty parameter to " << this->penalty_parameter << '\n';
      direction = this->solve_l1_relaxed_problem(statistics, current_iterate, this->penalty_parameter, warmstart_information);

      // recompute the linearized residual
      linearized_residual = this->original_model.compute_linearized_constraint_violation(direction.primals, current_iterate.evaluations.constraints,
            current_iterate.evaluations.constraint_jacobian, direction.primal_dual_step_length, Norm::L1);
      DEBUG << "Linearized infeasibility mk(dk): " << linearized_residual << "\n\n";
   }
   DEBUG << "Condition enforce_linearized_residual_sufficient_decrease is true\n";
   return direction;
}

bool l1Relaxation::linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual,
      double residual_lowest_violation) const {
   if (residual_lowest_violation <= this->parameters.residual_small_threshold) {
      return (linearized_residual <= this->parameters.residual_small_threshold);
   }
   const double linearized_residual_reduction = current_iterate.progress.infeasibility - linearized_residual;
   const double lowest_linearized_residual_reduction = current_iterate.progress.infeasibility - residual_lowest_violation;
   return (linearized_residual_reduction >= this->parameters.epsilon1 * lowest_linearized_residual_reduction);
}

Direction l1Relaxation::enforce_descent_direction_for_l1_merit(Statistics& statistics, Iterate& current_iterate, Direction& direction,
      const Direction& direction_lowest_violation, WarmstartInformation& warmstart_information) {
   while (0. < this->penalty_parameter && not this->is_descent_direction_for_l1_merit_function(current_iterate, direction, direction_lowest_violation)) {
      // decrease the penalty parameter and re-solve the problem
      this->penalty_parameter /= this->parameters.decrease_factor;
      DEBUG << "Further decrease the penalty parameter to " << this->penalty_parameter << '\n';
      direction = this->solve_l1_relaxed_problem(statistics, current_iterate, this->penalty_parameter, warmstart_information);
   }
   DEBUG << "Condition enforce_descent_direction_for_l1_merit is true\n";
   return direction;
}

bool l1Relaxation::is_descent_direction_for_l1_merit_function(const Iterate& current_iterate, const Direction& direction,
      const Direction& direction_lowest_violation) const {
   const double predicted_l1_merit_reduction = current_iterate.residuals.infeasibility - direction.subproblem_objective;
   const double lowest_decrease_objective = current_iterate.residuals.infeasibility - direction_lowest_violation.subproblem_objective;
   return (predicted_l1_merit_reduction >= this->parameters.epsilon2 * lowest_decrease_objective);
}

void l1Relaxation::compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& /*direction*/, double /*step_length*/) {
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed\n";
      this->globalization_strategy->reset();
      this->subproblem->subproblem_definition_changed = false;
   }
   // compute the progress measures for the current and trial iterates
   this->set_progress_measures(current_iterate);
   this->set_progress_measures(trial_iterate);

   trial_iterate.multipliers.objective = this->l1_relaxed_problem.get_objective_multiplier();
}

bool l1Relaxation::is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      double step_length) {
   this->subproblem->postprocess_iterate(this->l1_relaxed_problem, trial_iterate);
   this->compute_progress_measures(current_iterate, trial_iterate, direction, step_length);

   bool accept_iterate = false;
   if (direction.norm == 0.) {
      DEBUG << "Zero step acceptable\n\n";
      trial_iterate.evaluate_objective(this->original_model);
      accept_iterate = true;
   }
   else {
      // evaluate the predicted reduction
      ProgressMeasures predicted_reduction = this->compute_predicted_reduction_models(current_iterate, direction, step_length);

      // invoke the globalization strategy for acceptance
      accept_iterate = this->globalization_strategy->is_iterate_acceptable(statistics, trial_iterate, current_iterate.progress, trial_iterate.progress,
            predicted_reduction, this->penalty_parameter);
   }

   if (accept_iterate) {
      this->compute_primal_dual_residuals(this->original_model, this->feasibility_problem, trial_iterate);
      this->add_statistics(statistics, trial_iterate);
      this->check_exact_relaxation(trial_iterate);
   }
   return accept_iterate;
}

void l1Relaxation::set_progress_measures(Iterate& iterate) const {
   this->l1_relaxed_problem.set_infeasibility_measure(iterate, Norm::L1);
   this->l1_relaxed_problem.set_optimality_measure(iterate);
   this->subproblem->set_auxiliary_measure(this->l1_relaxed_problem, iterate);
}

ProgressMeasures l1Relaxation::compute_predicted_reduction_models(Iterate& current_iterate, const Direction& direction, double step_length) {
   return {
      this->l1_relaxed_problem.compute_predicted_infeasibility_reduction_model(current_iterate, direction, step_length, Norm::L1),
      this->subproblem->compute_predicted_optimality_reduction_model(this->l1_relaxed_problem, current_iterate, direction, step_length),
      this->subproblem->compute_predicted_auxiliary_reduction_model(this->l1_relaxed_problem, current_iterate, direction, step_length)
   };
}

double l1Relaxation::compute_complementarity_error(const std::vector<double>& primals, const std::vector<double>& constraints,
      const Multipliers& multipliers) const {
   return this->l1_relaxed_problem.compute_complementarity_error(primals, constraints, multipliers, Norm::L1);
}

void l1Relaxation::set_trust_region_radius(double trust_region_radius) {
   this->subproblem->set_trust_region_radius(trust_region_radius);
}

// for information purposes, check that l1 is an exact relaxation
void l1Relaxation::check_exact_relaxation(Iterate& iterate) const {
   const double norm_inf_multipliers = norm_inf(iterate.multipliers.constraints);
   if (0. < norm_inf_multipliers && this->penalty_parameter <= 1./norm_inf_multipliers) {
      DEBUG << "The value of the penalty parameter is consistent with an exact relaxation\n\n";
   }
}

void l1Relaxation::add_statistics(Statistics& statistics, const Iterate& trial_iterate) const {
   statistics.add_statistic("complementarity", trial_iterate.residuals.optimality_complementarity);
   statistics.add_statistic("stationarity", trial_iterate.residuals.optimality_stationarity);
   statistics.add_statistic("penalty param.", this->penalty_parameter);
   if (this->original_model.is_constrained()) {
      statistics.add_statistic("primal infeas.", trial_iterate.progress.infeasibility);
   }
}

size_t l1Relaxation::get_hessian_evaluation_count() const {
   return this->subproblem->get_hessian_evaluation_count();
}

size_t l1Relaxation::get_number_subproblems_solved() const {
   return this->subproblem->number_subproblems_solved;
}
