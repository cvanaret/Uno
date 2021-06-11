#include "l1Relaxation.hpp"
#include "GlobalizationStrategyFactory.hpp"

l1Relaxation::l1Relaxation(Problem& problem, Subproblem& subproblem, const std::map<std::string, std::string>& options) :
   ConstraintRelaxationStrategy(subproblem),
   globalization_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)),
   penalty_parameter(1.),
   parameters({10., 0.1, 0.1}) {
   // generate elastic variables to relax the constraints
   ActiveSetMethod::generate_elastic_variables_(problem, this->elastic_variables);

   // add nonnegative elastic variables to the subproblem
   this->subproblem.variables_bounds.resize(problem.number_variables + this->elastic_variables.size());
   for (size_t i = problem.number_variables; i < this->subproblem.variables_bounds.size(); i++) {
      this->subproblem.variables_bounds[i] = {0., INFINITY};
   }
}

Iterate l1Relaxation::initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
   /* initialize the subproblem */
   Iterate first_iterate = this->subproblem.evaluate_initial_point(problem, x, multipliers);
   this->globalization_strategy->initialize(statistics, first_iterate);

   // preallocate trial_iterate
   this->trial_primals_.resize(first_iterate.x.size());

   return first_iterate;
}

Direction l1Relaxation::compute_feasible_direction(const Problem& problem, Iterate& current_iterate, double trust_region_radius) {
   /*
   // preprocess the subproblem: scale the objective gradient and introduce the elastic variables
   this->preprocess_subproblem();

   // use Byrd's steering rules to update the penalty parameter and compute descent directions
   std::vector<Direction> directions = this->compute_byrd_steering_rule(problem, current_iterate, trust_region_radius);

   // postprocess the subproblem: remove the elastic variables
   for (Direction& direction: directions) {
      this->postprocess_direction(problem, direction);
   }
   return directions;
   */
   return this->subproblem.compute_direction(problem, current_iterate, trust_region_radius);
}

std::optional<Iterate> l1Relaxation::check_acceptance(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Direction& direction,
      double step_length) {
   // compute the trial iterate TODO do not reevaluate if ||d|| = 0
   add_vectors(current_iterate.x, direction.x, step_length, this->trial_primals_);
   Iterate trial_iterate(this->trial_primals_, direction.multipliers);
   double step_norm = step_length * direction.norm;
   this->subproblem.compute_optimality_measures(problem, trial_iterate);

   // check if subproblem definition changed
   if (this->subproblem.subproblem_definition_changed) {
      this->globalization_strategy->reset();
      this->subproblem.subproblem_definition_changed = false;
      this->subproblem.compute_optimality_measures(problem, current_iterate);
   }

   bool accept = false;
   if (step_norm == 0.) {
      accept = true;
   }
   else {
      // compute the predicted reduction (a mixture of the subproblem's and of the l1 relaxation's)
      double predicted_reduction = this->compute_predicted_reduction(problem, current_iterate, direction, step_length);
      accept = this->globalization_strategy->check_acceptance(statistics, current_iterate.progress, trial_iterate.progress, direction,
            predicted_reduction);
   }

   return (accept ? std::make_optional(trial_iterate) : std::nullopt);
}

void l1Relaxation::preprocess_subproblem() {
   // scale the objective gradient
   if (this->penalty_parameter == 0.) {
      this->subproblem.objective_gradient.clear();
   }
   else if (this->penalty_parameter < 1.) {
      for (auto& derivative : this->subproblem.objective_gradient) {
         derivative.second *= this->penalty_parameter;
      }
   }
   // add the positive elastic variables
   for (const auto& [j, i]: elastic_variables.positive) {
      this->subproblem.objective_gradient[i] = 1.;
      this->subproblem.constraints_jacobian[j][i] = -1.;
   }
   // add the negative elastic variables
   for (const auto& [j, i]: elastic_variables.negative) {
      this->subproblem.objective_gradient[i] = 1.;
      this->subproblem.constraints_jacobian[j][i] = 1.;
   }
}

std::vector<Direction> l1Relaxation::compute_byrd_steering_rule(const Problem& problem, Iterate& current_iterate, double trust_region_radius) {
   DEBUG << "penalty parameter: " << this->penalty_parameter << "\n";
   // TODO: pass penalty parameter to the Hessian and multipliers to the subproblem

   /* stage a: compute the step within trust region */
   // std::vector<Direction> directions = this->subproblem.compute_directions(problem, current_iterate, trust_region_radius);
   this->subproblem.generate(problem, current_iterate, this->penalty_parameter, trust_region_radius);
   Direction direction = this->subproblem.compute_direction(problem, current_iterate, trust_region_radius);

   /* penalty update: if penalty parameter is already 0, no need to decrease it */
   if (0. < this->penalty_parameter) {
      /* check infeasibility */
      double linearized_residual = this->compute_linearized_constraint_residual(direction.x);
      DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";

      // if problem had to be relaxed
      if (linearized_residual != 0.) {
         double current_penalty_parameter = this->penalty_parameter;

         /* stage c: solve the ideal l1 penalty problem with a zero penalty (no objective) */
         DEBUG << "Compute ideal solution:\n";
         this->subproblem.generate(problem, current_iterate, 0., trust_region_radius);
         Direction ideal_direction = this->subproblem.compute_direction(problem, current_iterate, trust_region_radius);

         /* compute the ideal error (with a zero penalty parameter) */
         double ideal_error = this->compute_error(problem, current_iterate, ideal_direction.multipliers, 0.);
         DEBUG << "Ideal error: " << ideal_error << "\n";

         if (ideal_error == 0.) {
            /* stage f: update the penalty parameter */
            this->penalty_parameter = 0.;
         }
         else {
            double ideal_linearized_residual = this->compute_linearized_constraint_residual(ideal_direction.x);
            DEBUG << "Linearized residual mk(dk): " << ideal_linearized_residual << "\n\n";

            /* decrease penalty parameter to satisfy 2 conditions */
            bool condition1 = false, condition2 = false;
            while (!condition2) {
               if (!condition1) {
                  /* stage d: reach a fraction of the ideal decrease */
                  if ((ideal_linearized_residual == 0. && linearized_residual == 0) || (ideal_linearized_residual != 0. &&
                     current_iterate.progress.feasibility - linearized_residual >= this->parameters.epsilon1 *
                     (current_iterate.progress.feasibility - ideal_linearized_residual))) {
                     condition1 = true;
                     DEBUG << "Condition 1 is true\n";
                  }
               }
               /* stage e: further decrease penalty parameter if necessary */
               if (condition1 && current_iterate.progress.feasibility - direction.objective >=
                                 this->parameters.epsilon2 * (current_iterate.progress.feasibility - ideal_direction.objective)) {
                  condition2 = true;
                  DEBUG << "Condition 2 is true\n";
               }
               if (!condition2) {
                  this->penalty_parameter /= this->parameters.tau;
                  if (this->penalty_parameter < 1e-10) {
                     this->penalty_parameter = 0.;
                     condition2 = true;
                  }
                  else {
                     DEBUG << "\nAttempting to solve with penalty parameter " << this->penalty_parameter << "\n";
                     this->subproblem.generate(problem, current_iterate, this->penalty_parameter, trust_region_radius);
                     direction = this->subproblem.compute_direction(problem, current_iterate, trust_region_radius);

                     linearized_residual = this->compute_linearized_constraint_residual(direction.x);
                     DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";
                  }
               }
            }

            /* stage f: update the penalty parameter */
            double updated_penalty_parameter = this->penalty_parameter;
            double term = ideal_error / std::max(1., current_iterate.progress.feasibility);
            this->penalty_parameter = std::min(this->penalty_parameter, term * term);
            if (this->penalty_parameter < updated_penalty_parameter) {
               this->subproblem.generate(problem, current_iterate, this->penalty_parameter, trust_region_radius);
               direction = this->subproblem.compute_direction(problem, current_iterate, trust_region_radius);
            }

         }

         if (this->penalty_parameter < current_penalty_parameter) {
            DEBUG << "\n*** Penalty parameter updated to " << this->penalty_parameter << "\n";
            // TODO this->subproblem_definition_changed = true;
            if (this->penalty_parameter == 0.) {
               direction = ideal_direction;
            }
         }
      }
   }
   direction.objective_multiplier = penalty_parameter;
   return std::vector<Direction>{direction};
}

double l1Relaxation::compute_linearized_constraint_residual(std::vector<double>& direction) {
   double residual = 0.;
   // l1 residual of the linearized constraints: sum of elastic variables
   for (std::pair<const int, int>& element: this->elastic_variables.positive) {
      int i = element.second;
      residual += direction[i];
   }
   for (std::pair<const int, int>& element: this->elastic_variables.negative) {
      int i = element.second;
      residual += direction[i];
   }
   return residual;
}

double l1Relaxation::compute_error(const Problem& problem, Iterate& iterate, Multipliers& multipliers, double penalty_parameter) {
   /* measure that combines KKT error and complementarity error */
   double error = 0.;

   /* KKT error */
   std::vector<double> lagrangian_gradient = iterate.lagrangian_gradient(problem, penalty_parameter, multipliers);
   error += norm_1(lagrangian_gradient);
   /* complementarity error */
   error += this->compute_complementarity_error(problem, iterate, multipliers);
   return error;
}

void l1Relaxation::postprocess_direction(const Problem& problem, Direction& direction) {
   /* remove p and n */
   direction.x.resize(problem.number_variables);
   direction.multipliers.lower_bounds.resize(problem.number_variables);
   direction.multipliers.upper_bounds.resize(problem.number_variables);
   direction.norm = norm_inf(direction.x);

   /* remove contribution of positive part variables */
   for (auto& [j, i]: elastic_variables.positive) {
      this->subproblem.objective_gradient.erase(i);
      this->subproblem.constraints_jacobian[j].erase(i);
   }
   /* remove contribution of negative part variables */
   for (auto& [j, i]: elastic_variables.negative) {
      this->subproblem.objective_gradient.erase(i);
      this->subproblem.constraints_jacobian[j].erase(i);
   }
}

double l1Relaxation::compute_predicted_reduction(const Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) {
   // compute the predicted reduction of the l1 relaxation as a postprocessing of the predicted reduction of the subproblem
   if (step_length == 1.) {
      return current_iterate.progress.feasibility + direction.predicted_reduction(step_length);
   }
   else {
      // determine the linearized constraint violation term: c(x_k) + alpha*\nabla c(x_k)^T d
      // TODO: creating a vector is not necessary!
      //auto linearized_constraint = [&](size_t j) {
      //   return current_iterate.constraints[j] + step_length * dot(direction.x, current_iterate.constraints_jacobian[j]);
      //}
      std::vector<double> linearized_constraints(current_iterate.constraints);
      for (size_t j = 0; j < current_iterate.constraints.size(); j++) {
         linearized_constraints[j] += step_length * dot(direction.x, current_iterate.constraints_jacobian[j]);
      }
      double linearized_constraint_violation = problem.compute_constraint_residual(linearized_constraints, L1_NORM);
      return current_iterate.progress.feasibility - linearized_constraint_violation + direction.predicted_reduction(step_length);
   }
}

/* complementary slackness error. Use abs/1e-8 to safeguard */
double l1Relaxation::compute_complementarity_error(const Problem& problem, Iterate& iterate, const Multipliers& multipliers) const {
   double error = 0.;
   /* bound constraints */
   for (size_t i = 0; i < problem.number_variables; i++) {
      if (-INFINITY < problem.variables_bounds[i].lb) {
         error += std::abs(iterate.multipliers.lower_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].lb));
      }
      if (problem.variables_bounds[i].ub < INFINITY) {
         error += std::abs(iterate.multipliers.upper_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].ub));
      }
   }
   /* general constraints */
   for (size_t j = 0; j < problem.number_constraints; j++) {
      double multiplier_j = multipliers.constraints[j];
      if (iterate.constraints[j] < problem.constraint_bounds[j].lb) {
         // violated lower: the multiplier is 1 at optimum
         error += std::abs((1. - multiplier_j) * (problem.constraint_bounds[j].lb - iterate.constraints[j]));
      }
      else if (problem.constraint_bounds[j].ub < iterate.constraints[j]) {
         // violated upper: the multiplier is -1 at optimum
         error += std::abs((1. + multiplier_j) * (iterate.constraints[j] - problem.constraint_bounds[j].ub));
      }
      else if (-INFINITY < problem.constraint_bounds[j].lb && 0. < multiplier_j) {
         error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraint_bounds[j].lb));
      }
      else if (problem.constraint_bounds[j].ub < INFINITY && multiplier_j < 0.) {
         error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraint_bounds[j].ub));
      }
   }
   return error;
}