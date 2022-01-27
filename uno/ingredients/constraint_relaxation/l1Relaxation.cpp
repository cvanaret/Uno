#include <cassert>
#include "l1Relaxation.hpp"
#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

l1Relaxation::l1Relaxation(const Problem& problem, const Options& options) :
      ConstraintRelaxationStrategy(problem, stod(options.at("l1_relaxation_initial_parameter")), options),
      globalization_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)),
      penalty_parameter(stod(options.at("l1_relaxation_initial_parameter"))),
      parameters({options.at("l1_relaxation_fixed_parameter") == "yes",
                  stod(options.at("l1_relaxation_decrease_factor")),
                  stod(options.at("l1_relaxation_epsilon1")),
                  stod(options.at("l1_relaxation_epsilon2")),
                  stod(options.at("l1_relaxation_small_threshold"))}),
      constraint_multipliers(problem.number_constraints) {
}

void l1Relaxation::initialize(Statistics& statistics, Iterate& first_iterate) {
   statistics.add_column("penalty param.", Statistics::double_width, 4);

   this->relaxed_problem.reset_elastic_variables(first_iterate);
   // initialize the subproblem
   this->subproblem->initialize(statistics, this->relaxed_problem, first_iterate);

   // compute the progress measures and the residuals of the initial point
   this->subproblem->compute_progress_measures(this->relaxed_problem, first_iterate);
   this->subproblem->compute_residuals(this->relaxed_problem, first_iterate, this->penalty_parameter);

   // initialize the globalization strategy
   this->globalization_strategy->initialize(statistics, first_iterate);
}

void l1Relaxation::create_current_subproblem(Iterate& current_iterate, double trust_region_radius) {
   // reset the elastic variables
   this->relaxed_problem.reset_elastic_variables(current_iterate);
   // scale the derivatives and introduce the elastic variables
   this->subproblem->build_current_subproblem(this->relaxed_problem, current_iterate, this->penalty_parameter, trust_region_radius);
   // set the multipliers of the violated constraints
   //this->set_multipliers(current_iterate, current_iterate.multipliers.constraints);
}

void l1Relaxation::set_multipliers(const Iterate& current_iterate, std::vector<double>& current_constraint_multipliers) {
   // the values {1, -1} are derived from the KKT conditions of the l1 problem
   for (size_t j = 0; j < this->original_problem.number_constraints; j++) {
      if (current_iterate.constraints[j] < this->original_problem.get_constraint_lower_bound(j)) { // lower bound infeasible
         current_constraint_multipliers[j] = 1.;
      }
      else if (this->original_problem.get_constraint_upper_bound(j) < current_iterate.constraints[j]) { // upper bound infeasible
         current_constraint_multipliers[j] = -1.;
      }
   }
   // otherwise, leave the multiplier as it is
}

Direction l1Relaxation::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) {
   DEBUG << "penalty parameter: " << this->penalty_parameter << "\n";
   // use Byrd's steering rules to update the penalty parameter and compute descent directions
   Direction direction = this->solve_with_steering_rule(statistics, current_iterate);
   return direction;
}

Direction l1Relaxation::solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
      const std::optional<std::vector<double>>& /*optional_phase_2_solution*/) {
   assert(0. < this->penalty_parameter && "l1Relaxation: the penalty parameter is already 0");
   Direction direction = this->resolve_subproblem(statistics, current_iterate, 0.);
   return direction;
}

Direction l1Relaxation::solve_subproblem(Statistics& statistics, Iterate& current_iterate, double current_penalty_parameter) {
   // solve the subproblem
   Direction direction = this->subproblem->solve(statistics, this->relaxed_problem, current_iterate);
   direction.objective_multiplier = current_penalty_parameter;
   DEBUG << "\n" << direction << "\n";
   assert(direction.status == OPTIMAL && "The subproblem was not solved to optimality");
   // enforce feasibility (by construction)
   if (direction.constraint_partition.has_value()) {
      const ConstraintPartition& constraint_partition = direction.constraint_partition.value();
      assert(constraint_partition.infeasible.empty() && "solve_subproblem: infeasible constraints found, although direction is feasible");
   }
   return direction;
}

Direction l1Relaxation::resolve_subproblem(Statistics& statistics, Iterate& current_iterate, double current_penalty_parameter) {
   this->relaxed_problem.set_objective_multiplier(current_penalty_parameter);
   // recompute the objective model with the current objective multiplier
   current_iterate.is_objective_gradient_computed = false;
   this->subproblem->build_objective_model(this->relaxed_problem, current_iterate, current_penalty_parameter);
   // solve the subproblem
   return this->solve_subproblem(statistics, current_iterate, current_penalty_parameter);
}

// Infeasibility detection and SQP methods for nonlinear optimization
// Richard H. Byrd, Frank E. Curtis and Jorge Nocedal
// https://epubs.siam.org/doi/pdf/10.1137/080738222
Direction l1Relaxation::solve_with_steering_rule(Statistics& statistics, Iterate& current_iterate) {
   // stage a: compute the step within trust region
   Direction direction = this->solve_subproblem(statistics, current_iterate, this->penalty_parameter);

   // penalty update: if penalty parameter is already 0 or fixed by the user, no need to decrease it
   if (0. < this->penalty_parameter && !this->parameters.fixed_parameter) {
      // check infeasibility
      double linearized_residual = this->relaxed_problem.compute_elastic_residual(current_iterate.x, direction.x);
      DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";

      // if the current direction is already feasible, terminate
      if (0. < linearized_residual) {
         const double current_penalty_parameter = this->penalty_parameter;

         // stage c: compute the lowest possible constraint violation (penalty parameter = 0)
         DEBUG << "Compute ideal solution (penalty parameter = 0):\n";
         Direction direction_lowest_violation = this->resolve_subproblem(statistics, current_iterate, 0.);
         const double residual_lowest_violation = this->relaxed_problem.compute_elastic_residual(current_iterate.x, direction_lowest_violation.x);
         DEBUG << "Lowest linearized residual mk(dk): " << residual_lowest_violation << "\n\n";

         // stage f: update the penalty parameter
         this->decrease_parameter_aggressively(current_iterate, direction_lowest_violation);
         DEBUG << "Penalty parameter aggressively set to " << this->penalty_parameter << "\n";
         if (this->penalty_parameter == 0.) {
            direction = direction_lowest_violation;
            linearized_residual = residual_lowest_violation;
         }
         else {
            if (this->penalty_parameter < current_penalty_parameter) {
               DEBUG << "Resolving the problem with penalty parameter " << this->penalty_parameter << "\n";
               direction = this->resolve_subproblem(statistics, current_iterate, this->penalty_parameter);
               linearized_residual = this->relaxed_problem.compute_elastic_residual(current_iterate.x, direction.x);
            }

            // further decrease penalty parameter to satisfy 2 conditions
            DEBUG << "Further decrease the penalty parameter\n";
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
                  this->penalty_parameter = std::max(0., this->penalty_parameter - this->parameters.small_threshold);
                  if (this->penalty_parameter == 0.) {
                     direction = direction_lowest_violation;
                     linearized_residual = residual_lowest_violation;
                     condition2 = true;
                  }
                  else {
                     DEBUG << "\nAttempting to solve with penalty parameter " << this->penalty_parameter << "\n";
                     direction = this->resolve_subproblem(statistics, current_iterate, this->penalty_parameter);
                     linearized_residual = this->relaxed_problem.compute_elastic_residual(current_iterate.x, direction.x);
                     DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";
                  }
               }
            }
         }

         if (this->penalty_parameter < current_penalty_parameter) {
            DEBUG << "\n*** Penalty parameter updated to " << this->penalty_parameter << "\n";
         }
      }
   }
   return direction;
}

bool l1Relaxation::linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual, double residual_lowest_violation) const {
   if (residual_lowest_violation == 0.) {
      return (linearized_residual == 0.);
   }
   const double linearized_residual_reduction = current_iterate.nonlinear_errors.constraints - linearized_residual;
   const double lowest_linearized_residual_reduction = current_iterate.nonlinear_errors.constraints - residual_lowest_violation;
   return (linearized_residual_reduction >= this->parameters.epsilon1 * lowest_linearized_residual_reduction);
}

void l1Relaxation::decrease_parameter_aggressively(Iterate& current_iterate, const Direction& direction_lowest_violation) {
   // compute the ideal error (with a zero penalty parameter)
   const double error_lowest_violation = l1Relaxation::compute_error(current_iterate, direction_lowest_violation.multipliers, 0.);
   DEBUG << "Ideal error: " << error_lowest_violation << "\n";

   const double scaled_error = error_lowest_violation / std::max(1., current_iterate.nonlinear_errors.constraints);
   const double scaled_error_square = scaled_error*scaled_error;
   this->penalty_parameter = std::min(this->penalty_parameter, scaled_error_square);
   this->penalty_parameter = std::max(0., this->penalty_parameter - this->parameters.small_threshold);
}

bool l1Relaxation::objective_sufficient_decrease(const Iterate& current_iterate, const Direction& direction,
      const Direction& direction_lowest_violation) const {
   const double decrease_objective = current_iterate.nonlinear_errors.constraints - direction.objective;
   const double lowest_decrease_objective = current_iterate.nonlinear_errors.constraints - direction_lowest_violation.objective;
   std::cout << decrease_objective << " >= " << this->parameters.epsilon2 << "*" << lowest_decrease_objective << " ?\n";
   return (decrease_objective >= this->parameters.epsilon2 * lowest_decrease_objective);
}

bool l1Relaxation::is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      PredictedReductionModel& predicted_reduction_model, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem->subproblem_definition_changed) {
      this->globalization_strategy->reset();
      this->subproblem->subproblem_definition_changed = false;
      this->subproblem->compute_progress_measures(this->original_problem, current_iterate);
   }

   bool accept = false;
   if (ConstraintRelaxationStrategy::is_small_step(direction)) {
      this->subproblem->compute_progress_measures(this->original_problem, trial_iterate);
      accept = true;
   }
   else {
      // evaluate the predicted reduction
      const double predicted_reduction = l1Relaxation::compute_predicted_reduction(this->original_problem, current_iterate, direction,
            predicted_reduction_model, step_length);
      // invoke the globalization strategy for acceptance
      trial_iterate.evaluate_constraints(this->original_problem);
      this->subproblem->compute_progress_measures(this->original_problem, trial_iterate);
      accept = this->globalization_strategy->check_acceptance(statistics, current_iterate.progress, trial_iterate.progress,
            this->penalty_parameter, predicted_reduction);
   }
   if (accept) {
      statistics.add_statistic("penalty param.", this->penalty_parameter);
      this->subproblem->compute_residuals(this->relaxed_problem, trial_iterate, direction.objective_multiplier);
   }
   return accept;
}

double l1Relaxation::compute_predicted_reduction(const Problem& problem, Iterate& current_iterate,
      const Direction& direction, PredictedReductionModel& predicted_reduction_model, double step_length) {
   // compute the predicted reduction of the l1 relaxation as a postprocessing of the predicted reduction of the subproblem
   if (step_length == 1.) {
      return current_iterate.nonlinear_errors.constraints + predicted_reduction_model.evaluate(step_length);
   }
   else {
      // determine the linearized constraint violation term: c(x_k) + alpha*\nabla c(x_k)^T d
      const auto residual_function = [&](size_t j) {
         const double component_j = current_iterate.constraints[j] + step_length * dot(direction.x, current_iterate.constraint_jacobian[j]);
         return problem.compute_constraint_violation(component_j, j);
      };
      const double linearized_constraint_violation = norm_1(residual_function, problem.number_constraints);
      return current_iterate.nonlinear_errors.constraints - linearized_constraint_violation + predicted_reduction_model.evaluate(step_length);
   }
}

// measure that combines KKT error and complementarity error
double l1Relaxation::compute_error(Iterate& current_iterate, const Multipliers& multiplier_displacements, double current_penalty_parameter) {
   // assemble the trial constraints multipliers
   for (size_t j = 0; j < this->original_problem.number_constraints; j++) {
      this->constraint_multipliers[j] = current_iterate.multipliers.constraints[j] + multiplier_displacements.constraints[j];
   }

   // complementarity error
   double error = Subproblem::compute_complementarity_error(this->original_problem, current_iterate, this->constraint_multipliers,
         multiplier_displacements.lower_bounds, multiplier_displacements.upper_bounds);
   // KKT error
   current_iterate.evaluate_lagrangian_gradient(this->original_problem, current_penalty_parameter, this->constraint_multipliers,
         multiplier_displacements.lower_bounds, multiplier_displacements.upper_bounds);
   error += norm_1(current_iterate.lagrangian_gradient);
   return error;
}

Direction l1Relaxation::compute_second_order_correction(Iterate& trial_iterate) {
   return this->subproblem->compute_second_order_correction(this->relaxed_problem, trial_iterate);
}

PredictedReductionModel l1Relaxation::generate_predicted_reduction_model(const Iterate& current_iterate, const Direction& direction) const {
   // the predicted reduction should be that of the original problem. It will then be post-processed in compute_predicted_reduction()
   return this->subproblem->generate_predicted_reduction_model(this->original_problem, current_iterate, direction);
}