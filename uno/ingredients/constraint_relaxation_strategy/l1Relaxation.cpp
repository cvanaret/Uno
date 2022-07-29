// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "l1Relaxation.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"
#include "tools/Range.hpp"

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
      relaxed_problem(model, options.get_double("l1_relaxation_initial_parameter"), options.get_double("l1_constraint_violation_coefficient"),
            options.get_bool("l1_use_proximal_term")),
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
      constraint_multipliers(model.number_constraints),
      lower_bound_multipliers(model.number_variables),
      upper_bound_multipliers(model.number_variables),
      statistics_penalty_parameter_column_order(options.get_int("statistics_penalty_parameter_column_order")) {
}

void l1Relaxation::initialize(Statistics& statistics, Iterate& first_iterate) {
   statistics.add_column("penalty param.", Statistics::double_width, this->statistics_penalty_parameter_column_order);

   this->subproblem->set_elastic_variables(this->relaxed_problem, first_iterate);
   // initialize the subproblem
   this->subproblem->initialize(statistics, this->relaxed_problem, first_iterate);

   // compute the progress measures and the residuals of the initial point
   first_iterate.nonlinear_progress.infeasibility = this->compute_infeasibility_measure(first_iterate);
   first_iterate.nonlinear_progress.optimality = this->subproblem->compute_optimality_measure(this->optimality_problem, first_iterate);
   this->compute_nonlinear_residuals(this->relaxed_problem, first_iterate);

   // initialize the globalization strategy
   this->globalization_strategy->initialize(first_iterate);
}

void l1Relaxation::set_multipliers(const Iterate& current_iterate, std::vector<double>& current_constraint_multipliers) {
   // the values {1, -1} are derived from the KKT conditions of the l1 problem
   for (size_t j = 0; j < this->optimality_problem.number_constraints; j++) {
      if (current_iterate.original_evaluations.constraints[j] < this->optimality_problem.get_constraint_lower_bound(j)) { // lower bound infeasible
         current_constraint_multipliers[j] = 1.;
      }
      else if (this->optimality_problem.get_constraint_upper_bound(j) < current_iterate.original_evaluations.constraints[j]) { // upper bound infeasible
         current_constraint_multipliers[j] = -1.;
      }
   }
   // otherwise, leave the multiplier as it is
}

Direction l1Relaxation::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) {
   // set the elastic variables
   this->subproblem->set_elastic_variables(this->relaxed_problem, current_iterate);

   // set the proximal coefficient
   this->relaxed_problem.set_proximal_coefficient(std::sqrt(this->penalty_parameter));
   this->relaxed_problem.set_proximal_reference_point(current_iterate.primals);

   // set the multipliers of the violated constraints
   // this->set_multipliers(current_iterate, current_iterate.multipliers.constraints);
   DEBUG << "Current iterate\n" << current_iterate << '\n';

   // use Byrd's steering rules to update the penalty parameter and compute a descent direction
   return this->solve_with_steering_rule(statistics, current_iterate);
}

Direction l1Relaxation::solve_subproblem(Statistics& statistics, Iterate& current_iterate, double current_penalty_parameter) {
   // update the objective model using the current penalty parameter
   DEBUG << "penalty parameter: " << current_penalty_parameter << "\n\n";
   this->relaxed_problem.set_objective_multiplier(current_penalty_parameter);

   // solve the subproblem
   Direction direction = this->subproblem->solve(statistics, this->relaxed_problem, current_iterate);
   direction.objective_multiplier = current_penalty_parameter;
   direction.norm = norm_inf(direction.primals, Range(this->model.number_variables));
   DEBUG << direction << '\n';
   assert(direction.status == Status::OPTIMAL && "The subproblem was not solved to optimality");
   return direction;
}

Direction l1Relaxation::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate) {
   assert(0. < this->penalty_parameter && "l1Relaxation: the penalty parameter is already 0");
   Direction direction = this->solve_subproblem(statistics, current_iterate, 0.);
   return direction;
}

Direction l1Relaxation::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate, const std::vector<double>& initial_point) {
   this->subproblem->set_initial_point(initial_point);
   return this->solve_feasibility_problem(statistics, current_iterate);
}

Direction l1Relaxation::solve_with_steering_rule(Statistics& statistics, Iterate& current_iterate) {
   // stage a: compute the step within trust region
   Direction direction = this->solve_subproblem(statistics, current_iterate, this->penalty_parameter);

   // penalty update: if penalty parameter is already 0 or fixed by the user, no need to decrease it
   if (0. < this->penalty_parameter && !this->parameters.fixed_parameter) {
      // check infeasibility
      double linearized_residual = this->relaxed_problem.compute_linearized_constraint_violation(current_iterate.primals, direction.primals);
      DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";

      // if the current direction is already feasible, terminate
      if (0. < linearized_residual) {
         const double current_penalty_parameter = this->penalty_parameter;

         // stage c: compute the lowest possible constraint violation (penalty parameter = 0)
         DEBUG << "Compute ideal solution by solving the feasibility problem:\n";
         Direction direction_lowest_violation = this->solve_subproblem(statistics, current_iterate, 0.);
         const double residual_lowest_violation = this->relaxed_problem
               .compute_linearized_constraint_violation(current_iterate.primals, direction_lowest_violation.primals);
         DEBUG << "Lowest linearized residual mk(dk): " << residual_lowest_violation << '\n';

         // stage f: update the penalty parameter
         this->decrease_parameter_aggressively(current_iterate, direction_lowest_violation);
         if (this->penalty_parameter == 0.) {
            direction = direction_lowest_violation;
            linearized_residual = residual_lowest_violation;
         }
         else {
            if (this->penalty_parameter < current_penalty_parameter) {
               DEBUG << "Resolving the subproblem (penalty parameter aggressively reduced to " << this->penalty_parameter << ")\n";
               direction = this->solve_subproblem(statistics, current_iterate, this->penalty_parameter);
               linearized_residual = this->relaxed_problem.compute_linearized_constraint_violation(current_iterate.primals, direction.primals);
            }

            // further decrease penalty parameter to satisfy 2 conditions
            bool condition1 = false, condition2 = false;
            while (!condition2) {
               if (!condition1) {
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
               if (!condition2) {
                  this->penalty_parameter /= this->parameters.decrease_factor;
                  if (this->penalty_parameter < this->parameters.small_threshold) {
                     this->penalty_parameter = 0.;
                  }
                  DEBUG << "Further decrease the penalty parameter to " << this->penalty_parameter << '\n';
                  if (this->penalty_parameter == 0.) {
                     direction = direction_lowest_violation;
                     linearized_residual = residual_lowest_violation;
                     condition2 = true;
                  }
                  else {
                     DEBUG << "Resolving the subproblem (penalty parameter reduced to " << this->penalty_parameter << ")\n";
                     direction = this->solve_subproblem(statistics, current_iterate, this->penalty_parameter);
                     linearized_residual = this->relaxed_problem.compute_linearized_constraint_violation(current_iterate.primals, direction.primals);
                     DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";
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
   const double linearized_residual_reduction = current_iterate.constraint_violation - linearized_residual;
   const double lowest_linearized_residual_reduction = current_iterate.constraint_violation - residual_lowest_violation;
   return (linearized_residual_reduction >= this->parameters.epsilon1 * lowest_linearized_residual_reduction);
}

void l1Relaxation::decrease_parameter_aggressively(Iterate& current_iterate, const Direction& direction_lowest_violation) {
   // compute the ideal error (with a zero penalty parameter)
   const double error_lowest_violation = l1Relaxation::compute_error(current_iterate, direction_lowest_violation.multipliers);
   DEBUG << "Ideal error: " << error_lowest_violation << '\n';

   const double scaled_error = error_lowest_violation / std::max(1., current_iterate.constraint_violation);
   const double scaled_error_square = scaled_error*scaled_error;
   this->penalty_parameter = std::min(this->penalty_parameter, scaled_error_square);
   if (this->penalty_parameter < this->parameters.small_threshold) {
      this->penalty_parameter = 0.;
   }
}

bool l1Relaxation::objective_sufficient_decrease(const Iterate& current_iterate, const Direction& direction,
      const Direction& direction_lowest_violation) const {
   const double decrease_objective = current_iterate.constraint_violation - direction.objective;
   const double lowest_decrease_objective = current_iterate.constraint_violation - direction_lowest_violation.objective;
   return (decrease_objective >= this->parameters.epsilon2 * lowest_decrease_objective);
}

bool l1Relaxation::is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      PredictedOptimalityReductionModel& predicted_optimality_reduction_model, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed, the optimality measure is recomputed\n";
      current_iterate.nonlinear_progress.optimality = this->subproblem->compute_optimality_measure(this->optimality_problem, current_iterate);
      this->globalization_strategy->reset();
      this->subproblem->subproblem_definition_changed = false;
   }

   // compute the measures of progress for the trial iterate
   trial_iterate.nonlinear_progress.infeasibility = this->compute_infeasibility_measure(trial_iterate);
   trial_iterate.nonlinear_progress.optimality = this->subproblem->compute_optimality_measure(this->optimality_problem, trial_iterate);

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
      accept = this->globalization_strategy->is_acceptable(current_iterate.nonlinear_progress, trial_iterate.nonlinear_progress,
            this->penalty_parameter, predicted_reduction);
   }
   if (accept) {
      statistics.add_statistic("penalty param.", this->penalty_parameter);
      trial_iterate.evaluate_objective(this->model);
      this->compute_nonlinear_residuals(this->relaxed_problem, trial_iterate);
   }
   return accept;
}

PredictedOptimalityReductionModel l1Relaxation::generate_predicted_optimality_reduction_model(const Direction& direction) const {
   return this->subproblem->generate_predicted_optimality_reduction_model(this->relaxed_problem, direction);
}

// measure that combines KKT error and complementarity error
double l1Relaxation::compute_error(Iterate& current_iterate, const Multipliers& multiplier_displacements) {
   // assemble the trial multipliers
   for (size_t j = 0; j < this->optimality_problem.number_constraints; j++) {
      this->constraint_multipliers[j] = current_iterate.multipliers.constraints[j] + multiplier_displacements.constraints[j];
   }
   for (size_t i = 0; i < this->optimality_problem.number_variables; i++) {
      this->lower_bound_multipliers[i] = current_iterate.multipliers.lower_bounds[i] + multiplier_displacements.lower_bounds[i];
      this->upper_bound_multipliers[i] = current_iterate.multipliers.upper_bounds[i] + multiplier_displacements.upper_bounds[i];
   }

   // complementarity error
   current_iterate.evaluate_constraints(this->model);
   double error = this->model.compute_complementarity_error(current_iterate.primals, current_iterate.original_evaluations.constraints,
         this->constraint_multipliers, this->lower_bound_multipliers, this->upper_bound_multipliers);
   // KKT error
   current_iterate.evaluate_lagrangian_gradient(this->model, this->relaxed_problem.get_objective_multiplier(),
         this->constraint_multipliers, this->lower_bound_multipliers, this->upper_bound_multipliers);
   error += norm_1(current_iterate.lagrangian_gradient);
   return error;
}

void l1Relaxation::set_variable_bounds(const Iterate& current_iterate, double trust_region_radius) {
   this->subproblem->set_variable_bounds(this->relaxed_problem, current_iterate, trust_region_radius);
}

Direction l1Relaxation::compute_second_order_correction(Iterate& trial_iterate) {
   return this->subproblem->compute_second_order_correction(this->relaxed_problem, trial_iterate);
}

double l1Relaxation::compute_infeasibility_measure(Iterate& iterate) {
   iterate.evaluate_constraints(this->model);
   return this->model.compute_constraint_violation(iterate.original_evaluations.constraints, L1_NORM);
}

void l1Relaxation::register_accepted_iterate(Iterate& iterate) {
   // TODO check problem
   this->subproblem->postprocess_accepted_iterate(this->relaxed_problem, iterate);
   // check that l1 is an exact relaxation
   if (this->penalty_parameter <= 1./norm_inf(iterate.multipliers.constraints)) {
      DEBUG << "The value of the penalty parameter is consistent with an exact relaxation\n";
   }
}

size_t l1Relaxation::get_hessian_evaluation_count() const {
   return this->subproblem->get_hessian_evaluation_count();
}

size_t l1Relaxation::get_number_subproblems_solved() const {
   return this->subproblem->number_subproblems_solved;
}