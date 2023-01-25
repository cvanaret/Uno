// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "l1Relaxation.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

/*
 * Infeasibility detection and SQP methods for nonlinear optimization
 * Richard H. Byrd, Frank E. Curtis and Jorge Nocedal
 * http://epubs.siam.org/doi/pdf/10.1137/080738222
 */

l1Relaxation::l1Relaxation(const Model& model, const Options& options) :
      ConstraintRelaxationStrategy(model, options),
      // create the optimality problem
      optimality_problem(model),
      // create the relaxed problem by introducing elastic variables
      relaxed_problem(model, options.get_double("l1_relaxation_initial_parameter")),
      subproblem(SubproblemFactory::create(this->relaxed_problem.number_variables, this->relaxed_problem.number_constraints,
            this->relaxed_problem.get_maximum_number_hessian_nonzeros(), options)),
      globalization_strategy(GlobalizationStrategyFactory::create(options.get_string("strategy"), options)),
      penalty_parameter(options.get_double("l1_relaxation_initial_parameter")),
      parameters({
         options.get_bool("l1_relaxation_fixed_parameter"),
         options.get_double("l1_relaxation_decrease_factor"),
         options.get_double("l1_relaxation_epsilon1"),
         options.get_double("l1_relaxation_epsilon2"),
         options.get_double("l1_relaxation_small_threshold")
      }),
      trial_multipliers(this->relaxed_problem.number_variables, model.number_constraints),
      statistics_penalty_parameter_column_order(options.get_int("statistics_penalty_parameter_column_order")) {
}

void l1Relaxation::initialize(Statistics& statistics, Iterate& first_iterate) {
   statistics.add_column("penalty param.", Statistics::double_width, this->statistics_penalty_parameter_column_order);

   this->subproblem->set_elastic_variable_values(this->relaxed_problem, first_iterate);
   // initialize the subproblem
   this->subproblem->initialize(statistics, this->relaxed_problem, first_iterate);

   // compute the progress measures and the residuals of the initial point
   this->set_infeasibility_measure(first_iterate);
   this->set_scaled_optimality_measure(first_iterate);
   this->subproblem->set_unscaled_optimality_measure(this->relaxed_problem, first_iterate);

   //ConstraintRelaxationStrategy::evaluate_functions(this->relaxed_problem, first_iterate);
   this->compute_primal_dual_residuals(this->relaxed_problem, first_iterate);

   // initialize the globalization strategy
   this->globalization_strategy->initialize(first_iterate);
}

Direction l1Relaxation::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) {
   DEBUG << "Current iterate\n" << current_iterate << '\n';

   // use Byrd's steering rules to update the penalty parameter and compute a descent direction
   return this->solve_with_steering_rule(statistics, current_iterate);
}

Direction l1Relaxation::solve_subproblem(Statistics& statistics, Iterate& current_iterate, double current_penalty_parameter) {
   DEBUG << "Solving the subproblem with penalty parameter " << current_penalty_parameter << "\n\n";
   this->relaxed_problem.set_objective_multiplier(current_penalty_parameter);

   // solve the subproblem
   Direction direction = this->subproblem->solve(statistics, this->relaxed_problem, current_iterate);
   direction.objective_multiplier = current_penalty_parameter;
   direction.norm = norm_inf(direction.primals, Range(this->optimality_problem.number_variables));
   DEBUG << direction << '\n';
   assert(direction.status == SubproblemStatus::OPTIMAL && "The subproblem was not solved to optimality");
   return direction;
}

Direction l1Relaxation::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate) {
   assert(0. < this->penalty_parameter && "l1Relaxation: the penalty parameter is already 0");
   this->subproblem->initialize_feasibility_problem();
   return this->solve_subproblem(statistics, current_iterate, 0.);
}

Direction l1Relaxation::compute_second_order_correction(Iterate& trial_iterate) {
   // evaluate the constraints for the second-order correction
   Direction soc_direction = this->subproblem->compute_second_order_correction(this->relaxed_problem, trial_iterate);
   soc_direction.objective_multiplier = this->penalty_parameter;
   soc_direction.norm = norm_inf(soc_direction.primals, Range(this->optimality_problem.number_variables));
   DEBUG << soc_direction << '\n';
   return soc_direction;
}

Direction l1Relaxation::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate, const std::vector<double>& initial_point) {
   this->subproblem->set_initial_point(initial_point);
   return this->solve_feasibility_problem(statistics, current_iterate);
}

Direction l1Relaxation::solve_with_steering_rule(Statistics& statistics, Iterate& current_iterate) {
   // stage a: compute the step within trust region
   Direction direction = this->solve_subproblem(statistics, current_iterate, this->penalty_parameter);

   // penalty update: if penalty parameter is already 0 or fixed by the user, no need to decrease it
   if (0. < this->penalty_parameter && not this->parameters.fixed_parameter) {
      // check infeasibility
      double linearized_residual = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model, current_iterate,
            direction, 1.);
      DEBUG << "Linearized infeasibility mk(dk): " << linearized_residual << "\n\n";

      // if the current direction is already feasible, terminate
      if (0. < linearized_residual) {
         const double current_penalty_parameter = this->penalty_parameter;

         // stage c: compute the lowest possible constraint violation (penalty parameter = 0)
         DEBUG << "Compute ideal solution by solving the feasibility problem:\n";
         Direction direction_lowest_violation = this->solve_subproblem(statistics, current_iterate, 0.);
         const double residual_lowest_violation = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model,
               current_iterate, direction_lowest_violation, 1.);
         DEBUG << "Lowest linearized infeasibility mk(dk): " << residual_lowest_violation << '\n';

         // stage f: update the penalty parameter
         this->decrease_parameter_aggressively(current_iterate, direction_lowest_violation);
         DEBUG << "Further aggressively decrease the penalty parameter to " << this->penalty_parameter << '\n';
         if (this->penalty_parameter == 0.) {
            direction = direction_lowest_violation;
            linearized_residual = residual_lowest_violation;
         }
         else {
            if (this->penalty_parameter < current_penalty_parameter) {
               direction = this->solve_subproblem(statistics, current_iterate, this->penalty_parameter);
               linearized_residual = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model, current_iterate,
                     direction, 1.);
            }

            // further decrease penalty parameter to satisfy 2 conditions
            bool condition1 = false, condition2 = false;
            while (not condition2) {
               if (not condition1) {
                  // stage d: reach a fraction of the ideal decrease
                  if (this->linearized_residual_sufficient_decrease(current_iterate, linearized_residual, residual_lowest_violation)) {
                     condition1 = true;
                     DEBUG << "Condition 1 is true\n";
                  }
               }
               // stage e: further decrease penalty parameter if necessary
               if (condition1 && this->objective_sufficient_decrease(current_iterate, direction, direction_lowest_violation)) {
                  condition2 = true;
                  DEBUG << "Condition 2 is true\n";
               }
               if (not condition2) {
                  this->penalty_parameter /= this->parameters.decrease_factor;
                  DEBUG << "Further decrease the penalty parameter to " << this->penalty_parameter << '\n';
                  if (this->penalty_parameter == 0.) {
                     direction = direction_lowest_violation;
                     linearized_residual = residual_lowest_violation;
                     condition2 = true;
                  }
                  else {
                     direction = this->solve_subproblem(statistics, current_iterate, this->penalty_parameter);
                     linearized_residual = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model, current_iterate,
                           direction, 1.);
                     DEBUG << "Linearized infeasibility mk(dk): " << linearized_residual << "\n\n";
                  }
               }
            }
         }

         if (this->penalty_parameter < current_penalty_parameter) {
            DEBUG << "\nPenalty parameter updated to " << this->penalty_parameter << '\n';
         }
         DEBUG << '\n';
      }
   }
   return direction;
}

bool l1Relaxation::linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual, double residual_lowest_violation) const {
   if (residual_lowest_violation == 0.) {
      return (linearized_residual == 0.);
   }
   const double linearized_residual_reduction = current_iterate.residuals.infeasibility - linearized_residual;
   const double lowest_linearized_residual_reduction = current_iterate.residuals.infeasibility - residual_lowest_violation;
   return (linearized_residual_reduction >= this->parameters.epsilon1 * lowest_linearized_residual_reduction);
}

void l1Relaxation::decrease_parameter_aggressively(Iterate& current_iterate, const Direction& direction_lowest_violation) {
   // assemble the trial multipliers
   add_vectors(current_iterate.multipliers.constraints, direction_lowest_violation.multipliers.constraints, 1., this->trial_multipliers.constraints);
   add_vectors(current_iterate.multipliers.lower_bounds, direction_lowest_violation.multipliers.lower_bounds, 1., this->trial_multipliers.lower_bounds);
   add_vectors(current_iterate.multipliers.upper_bounds, direction_lowest_violation.multipliers.upper_bounds, 1., this->trial_multipliers.upper_bounds);

   // the ideal error (with penalty parameter = 0) must make sense: there must be at least a nonzero dual to avoid trivial FJ points
   if (this->trial_multipliers.not_all_zero(this->original_model.number_variables, this->small_duals_threshold)) {
      // compute the ideal error (with a zero penalty parameter)
      const double error_lowest_violation = l1Relaxation::compute_dual_error(current_iterate);
      DEBUG << "Ideal error: " << error_lowest_violation << '\n';
      const double scaled_error = error_lowest_violation / std::max(1., current_iterate.residuals.infeasibility);
      const double scaled_error_square = scaled_error * scaled_error;
      this->penalty_parameter = std::min(this->penalty_parameter, scaled_error_square);
   }
   else {
      WARNING << RED << "All multipliers are close to 0. The dual error shouldn't be used in l1Relaxation" << RESET << '\n';
   }
}

bool l1Relaxation::objective_sufficient_decrease(const Iterate& current_iterate, const Direction& direction,
      const Direction& direction_lowest_violation) const {
   const double decrease_objective = current_iterate.residuals.infeasibility - direction.subproblem_objective;
   const double lowest_decrease_objective = current_iterate.residuals.infeasibility - direction_lowest_violation.subproblem_objective;
   return (decrease_objective >= this->parameters.epsilon2 * lowest_decrease_objective);
}

void l1Relaxation::compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& /*direction*/) {
   // refresh the progress measures for the current iterate
   if (this->subproblem->unscaled_optimality_measure_changed) {
      DEBUG << "The subproblem definition changed, the unscaled optimality measure is recomputed\n";
      this->subproblem->set_unscaled_optimality_measure(this->relaxed_problem, current_iterate);
      this->globalization_strategy->reset();
      this->subproblem->unscaled_optimality_measure_changed = false;
   }
   this->set_infeasibility_measure(current_iterate);
   this->set_scaled_optimality_measure(current_iterate);

   // compute the progress measures for the trial iterate
   this->set_infeasibility_measure(trial_iterate);
   this->set_scaled_optimality_measure(trial_iterate);
   this->subproblem->set_unscaled_optimality_measure(this->relaxed_problem, trial_iterate);
}

bool l1Relaxation::is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      double step_length) {
   this->compute_progress_measures(current_iterate, trial_iterate, direction);

   bool accept = false;
   if (this->is_small_step(direction)) {
      DEBUG << "Small step acceptable\n";
      trial_iterate.evaluate_objective(this->original_model);
      accept = true;
   }
   else {
      // evaluate the predicted reduction
      PredictedReduction predicted_reduction = {
            this->generate_predicted_infeasibility_reduction_model(current_iterate, direction, step_length),
            this->generate_predicted_scaled_optimality_reduction_model(current_iterate, direction, step_length),
            this->subproblem->generate_predicted_unscaled_optimality_reduction_model(this->relaxed_problem, current_iterate, direction, step_length)
      };
      // invoke the globalization strategy for acceptance
      accept = this->globalization_strategy->is_iterate_acceptable(current_iterate.progress, trial_iterate.progress,
            predicted_reduction, this->penalty_parameter);
   }
   if (accept) {
      statistics.add_statistic("penalty param.", this->penalty_parameter);
      //ConstraintRelaxationStrategy::evaluate_functions(this->relaxed_problem, trial_iterate);
      this->compute_primal_dual_residuals(this->relaxed_problem, trial_iterate);
   }
   return accept;
}

void l1Relaxation::set_infeasibility_measure(Iterate& iterate) {
   if (0. < this->penalty_parameter) {
      // constraint violation
      iterate.evaluate_constraints(this->original_model);
      iterate.progress.infeasibility = this->original_model.compute_constraint_violation(iterate.evaluations.constraints, L1_NORM);
   }
   else {
      // 0
      iterate.progress.infeasibility = 0.;
   }
}

double l1Relaxation::generate_predicted_infeasibility_reduction_model(const Iterate& current_iterate, const Direction& direction, double step_length) const {
   if (0. < this->penalty_parameter) {
      const double current_constraint_violation = this->original_model.compute_constraint_violation(current_iterate.evaluations.constraints,
            Norm::L1_NORM);
      const double linearized_constraint_violation = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model,
            current_iterate, direction, step_length);
      return current_constraint_violation - linearized_constraint_violation;
      // "‖c(x)‖₁ - ‖c(x) + ∇c(x)^T (αd)‖₁"};
   }
   else {
      return 0.;
      // "0"};
   }
}

void l1Relaxation::set_scaled_optimality_measure(Iterate& iterate) {
   if (0. < this->penalty_parameter) {
      // scaled objective
      iterate.evaluate_objective(this->original_model);
      const double objective = iterate.evaluations.objective;
      iterate.progress.scaled_optimality = [=](double objective_multiplier) {
         return objective_multiplier*objective;
      };
   }
   else {
      // constraint violation
      iterate.evaluate_constraints(this->original_model);
      const double constraint_violation = this->original_model.compute_constraint_violation(iterate.evaluations.constraints, L1_NORM);
      iterate.progress.scaled_optimality = [=](double /*objective_multiplier*/) {
         return constraint_violation;
      };
   }
}

std::function<double (double)> l1Relaxation::generate_predicted_scaled_optimality_reduction_model(const Iterate& current_iterate,
      const Direction& direction, double step_length) const {
   if (0. < this->penalty_parameter) {
      // precompute expensive quantities
      const double directional_derivative = dot(direction.primals, current_iterate.evaluations.objective_gradient);
      return [=](double objective_multiplier) {
         return step_length * (-objective_multiplier*directional_derivative);
      };
      // "-ρ*∇f(x)^T (αd)"};
   }
   else {
      const double current_constraint_violation = this->original_model.compute_constraint_violation(current_iterate.evaluations.constraints,
            Norm::L1_NORM);
      const double linearized_constraint_violation = ConstraintRelaxationStrategy::compute_linearized_constraint_violation(this->original_model,
            current_iterate, direction, step_length);
      return [=](double /*objective_multiplier*/) {
         return current_constraint_violation - linearized_constraint_violation;
      };
      // "‖c(x)‖_1 - ‖c(x) + ∇c(x)^T (αd)‖_1"};
   }
}

// measure that combines KKT error and complementarity error
double l1Relaxation::compute_dual_error(Iterate& current_iterate) {
   // stationarity error
   ConstraintRelaxationStrategy::evaluate_lagrangian_gradient(this->original_model.number_variables, current_iterate, this->trial_multipliers,
         this->penalty_parameter);
   double error = norm_1(current_iterate.lagrangian_gradient.constraints_contribution);
   // complementarity error
   error += this->relaxed_problem.compute_feasibility_complementarity_error(this->original_model.number_variables, current_iterate.primals,
         current_iterate.evaluations.constraints, this->trial_multipliers);
   return error;
}

void l1Relaxation::set_trust_region_radius(double trust_region_radius) {
   this->subproblem->set_trust_region_radius(trust_region_radius);
}

void l1Relaxation::postprocess_accepted_iterate(Iterate& iterate) {
   this->subproblem->postprocess_accepted_iterate(this->relaxed_problem, iterate);

   // for information, check that l1 is an exact relaxation
   const double norm_inf_multipliers = norm_inf(iterate.multipliers.constraints);
   if (0. < norm_inf_multipliers && this->penalty_parameter <= 1./norm_inf_multipliers) {
      DEBUG << "The value of the penalty parameter is consistent with an exact relaxation\n\n";
   }
}

size_t l1Relaxation::get_hessian_evaluation_count() const {
   return this->subproblem->get_hessian_evaluation_count();
}

size_t l1Relaxation::get_number_subproblems_solved() const {
   return this->subproblem->number_subproblems_solved;
}
