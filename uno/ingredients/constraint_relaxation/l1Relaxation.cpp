#include <cassert>
#include "l1Relaxation.hpp"
#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

l1Relaxation::l1Relaxation(Problem& problem, const Options& options, bool use_trust_region) :
      ConstraintRelaxationStrategy(SubproblemFactory::create(problem, problem.number_variables + l1Relaxation::count_elastic_variables(problem),
            options.at("subproblem"), options, use_trust_region)),
      globalization_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)),
      penalty_parameter(stod(options.at("l1_relaxation_initial_parameter"))),
      parameters({10., 0.1, 0.1}) {

   // generate elastic variables to relax the constraints
   ConstraintRelaxationStrategy::generate_elastic_variables(problem, this->elastic_variables);
   // add bounds for nonnegative elastic variables to the subproblem
   this->set_elastic_bounds_in_subproblem(problem, this->elastic_variables.size());

   // elastic variables are temporary and are discarded when passed to mechanism
   this->number_variables -= l1Relaxation::count_elastic_variables(problem);
}

Iterate l1Relaxation::initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
   statistics.add_column("penalty param.", Statistics::double_width, 4);

   // initialize the subproblem
   Iterate first_iterate = this->subproblem->generate_initial_iterate(statistics, problem, x, multipliers);
   this->subproblem->compute_errors(problem, first_iterate, this->penalty_parameter);

   this->globalization_strategy->initialize(statistics, first_iterate);
   return first_iterate;
}

void l1Relaxation::generate_subproblem(const Problem& problem, Iterate& current_iterate, double trust_region_radius) {
   // preprocess the subproblem: scale the objective gradient and introduce the elastic variables
   // scale the objective gradient and (possibly) and Hessian
   this->subproblem->generate(problem, current_iterate, this->penalty_parameter, trust_region_radius);
   this->add_elastic_variables_to_subproblem(this->elastic_variables);
}

Direction l1Relaxation::compute_feasible_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   DEBUG << "penalty parameter: " << this->penalty_parameter << "\n";
   // use Byrd's steering rules to update the penalty parameter and compute descent directions
   Direction direction = this->compute_byrd_steering_rule(statistics, problem, current_iterate);

   // remove the temporary elastic variables
   this->remove_elastic_variables(problem, direction);
   return direction;
}

Direction l1Relaxation::solve_feasibility_problem(Statistics& statistics, const Problem& problem, Iterate& current_iterate,
      Direction& /*direction*/) {
   this->update_objective_multiplier(problem, current_iterate, 0.);
   return this->subproblem->solve(statistics, problem, current_iterate);
}

bool l1Relaxation::is_acceptable(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Iterate& trial_iterate,
      const Direction& direction, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem->subproblem_definition_changed) {
      this->globalization_strategy->reset();
      this->subproblem->subproblem_definition_changed = false;
      this->subproblem->compute_progress_measures(problem, current_iterate);
   }
   const double step_norm = step_length * direction.norm;

   bool accept = false;
   if (step_norm == 0.) {
      accept = true;
   }
   else {
      trial_iterate.compute_constraints(problem);
      // compute the predicted reduction (a mixture of the subproblem's and of the l1 relaxation's)
      const double predicted_reduction = this->compute_predicted_reduction(problem, current_iterate, direction, step_length);
      // invoke the globalization strategy for acceptance
      accept = this->globalization_strategy->check_acceptance(statistics, current_iterate.progress, trial_iterate.progress,
            this->penalty_parameter, predicted_reduction);
   }
   if (accept) {
      statistics.add_statistic("penalty param.", this->penalty_parameter);

      // compute the residuals
      trial_iterate.compute_objective(problem);
      this->subproblem->compute_errors(problem, trial_iterate, direction.objective_multiplier);
   }
   return accept;
}

void l1Relaxation::update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) {
   this->subproblem->update_objective_multiplier(problem, current_iterate, objective_multiplier);
   // add the positive elastic variables
   for (const auto& element: elastic_variables.positive) {
      const size_t i = element.second;
      this->subproblem->objective_gradient[i] = 1.;
   }
   // add the negative elastic variables
   for (const auto& element: elastic_variables.negative) {
      const size_t i = element.second;
      this->subproblem->objective_gradient[i] = 1.;
   }
}

Direction l1Relaxation::solve_subproblem(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   Direction direction = this->subproblem->solve(statistics, problem, current_iterate);
   assert(direction.constraint_partition.infeasible.empty() && "Infeasible constraints found, although direction is feasible");
   direction.objective_multiplier = this->penalty_parameter;
   return direction;
}

Direction l1Relaxation::solve_subproblem(Statistics& statistics, const Problem& problem, Iterate& current_iterate, double objective_multiplier) {
   this->update_objective_multiplier(problem, current_iterate, objective_multiplier);
   Direction direction = this->subproblem->solve(statistics, problem, current_iterate);
   assert(direction.constraint_partition.infeasible.empty() && "Infeasible constraints found, although direction is feasible");
   direction.objective_multiplier = objective_multiplier;
   return direction;
}

Direction l1Relaxation::compute_byrd_steering_rule(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   /* stage a: compute the step within trust region */
   Direction direction = this->solve_subproblem(statistics, problem, current_iterate);

   /* penalty update: if penalty parameter is already 0, no need to decrease it */
   if (0. < this->penalty_parameter) {
      /* check infeasibility */
      double linearized_residual = this->compute_linearized_constraint_residual(direction.x);
      DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";

      // if problem had to be relaxed
      if (linearized_residual != 0.) {
         const double current_penalty_parameter = this->penalty_parameter;

         /* stage c: solve the ideal l1 penalty problem with a zero penalty (no objective) */
         DEBUG << "Compute ideal solution (param = 0):\n";
         Direction ideal_direction = this->solve_subproblem(statistics, problem, current_iterate, 0.);

         /* compute the ideal error (with a zero penalty parameter) */
         const double ideal_error = l1Relaxation::compute_error(problem, current_iterate, ideal_direction.multipliers, 0.);
         DEBUG << "Ideal error: " << ideal_error << "\n";
         if (ideal_error == 0.) {
            /* stage f: update the penalty parameter */
            this->penalty_parameter = 0.;
         }
         else {
            const double ideal_linearized_residual = this->compute_linearized_constraint_residual(ideal_direction.x);
            DEBUG << "Linearized residual mk(dk): " << ideal_linearized_residual << "\n\n";

            /* decrease penalty parameter to satisfy 2 conditions */
            bool condition1 = false, condition2 = false;
            while (!condition2) {
               if (!condition1) {
                  /* stage d: reach a fraction of the ideal decrease */
                  if ((ideal_linearized_residual == 0. && linearized_residual == 0) || (ideal_linearized_residual != 0. &&
                                                                                        current_iterate.progress.infeasibility -
                                                                                        linearized_residual >= this->parameters.epsilon1 *
                                                                                                               (current_iterate.progress
                                                                                                                      .infeasibility -
                                                                                                                ideal_linearized_residual))) {
                     condition1 = true;
                     DEBUG << "Condition 1 is true\n";
                  }
               }
               /* stage e: further decrease penalty parameter if necessary */
               if (condition1 && current_iterate.progress.infeasibility - direction.objective >=
                                 this->parameters.epsilon2 * (current_iterate.progress.infeasibility - ideal_direction.objective)) {
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
                     direction = this->solve_subproblem(statistics, problem, current_iterate, this->penalty_parameter);

                     linearized_residual = this->compute_linearized_constraint_residual(direction.x);
                     DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";
                  }
               }
            }

            /* stage f: update the penalty parameter */
            double updated_penalty_parameter = this->penalty_parameter;
            double term = ideal_error / std::max(1., current_iterate.progress.infeasibility);
            this->penalty_parameter = std::min(this->penalty_parameter, term * term);
            if (this->penalty_parameter < updated_penalty_parameter) {
               direction = this->solve_subproblem(statistics, problem, current_iterate, this->penalty_parameter);
            }

         } // end else

         if (this->penalty_parameter < current_penalty_parameter) {
            DEBUG << "\n*** Penalty parameter updated to " << this->penalty_parameter << "\n";
            this->globalization_strategy->reset();
            if (this->penalty_parameter == 0.) {
               direction = ideal_direction;
            }
         }
      }
   }
   return direction;
}

size_t l1Relaxation::count_elastic_variables(const Problem& problem) {
   size_t number_variables = 0;
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (-INFINITY < problem.constraint_bounds[j].lb) {
         number_variables++;
      }
      if (problem.constraint_bounds[j].ub < INFINITY) {
         number_variables++;
      }
   }
   return number_variables;
}

double l1Relaxation::compute_linearized_constraint_residual(std::vector<double>& direction) {
   double residual = 0.;
   // l1 residual of the linearized constraints: sum of elastic variables
   for (const auto& element: this->elastic_variables.positive) {
      size_t i = element.second;
      residual += direction[i];
   }
   for (const auto& element: this->elastic_variables.negative) {
      size_t i = element.second;
      residual += direction[i];
   }
   return residual;
}

/* measure that combines KKT error and complementarity error */
double l1Relaxation::compute_error(const Problem& problem, Iterate& iterate, Multipliers& multipliers, double penalty_parameter) const {
   /* complementarity error */
   double error = this->subproblem->compute_complementarity_error(problem, iterate, multipliers);
   /* KKT error */
   std::vector<double> lagrangian_gradient = iterate.lagrangian_gradient(problem, penalty_parameter, multipliers);
   error += norm_1(lagrangian_gradient);
   return error;
}

void l1Relaxation::remove_elastic_variables(const Problem& problem, Direction& direction) {
   // compute set of satisfied/violated constraints

   /* remove p and n */
   // TODO change that!!!
   direction.x.resize(problem.number_variables);
   direction.multipliers.lower_bounds.resize(problem.number_variables);
   direction.multipliers.upper_bounds.resize(problem.number_variables);
   direction.norm = norm_inf(direction.x);

   /* remove contribution of positive part variables */
   for (auto&[j, i]: elastic_variables.positive) {
      this->subproblem->objective_gradient.erase(i);
      this->subproblem->constraints_jacobian[j].erase(i);
   }
   /* remove contribution of negative part variables */
   for (auto&[j, i]: elastic_variables.negative) {
      this->subproblem->objective_gradient.erase(i);
      this->subproblem->constraints_jacobian[j].erase(i);
   }
}

double l1Relaxation::compute_predicted_reduction(const Problem& problem, Iterate& current_iterate, const Direction& direction, double step_length) {
   // compute the predicted reduction of the l1 relaxation as a postprocessing of the predicted reduction of the subproblem
   if (step_length == 1.) {
      return current_iterate.progress.infeasibility + this->subproblem->compute_predicted_reduction(direction, step_length);;
   }
   else {
      // determine the linearized constraint violation term: c(x_k) + alpha*\nabla c(x_k)^T d
      auto residual_function = [&](size_t j) {
         const double component_j = current_iterate.constraints[j] + step_length * dot(direction.x, current_iterate.constraints_jacobian[j]);
         return problem.compute_constraint_violation(component_j, j);
      };
      const double linearized_constraint_violation = norm_1(residual_function, problem.number_constraints);
      return current_iterate.progress.infeasibility - linearized_constraint_violation + this->subproblem->compute_predicted_reduction(direction,
            step_length);
   }
}

void l1Relaxation::recover_l1qp_active_set_(const Problem& problem, const Direction& direction) {
   // remove extra variables p and n
   for (size_t i = problem.number_variables; i < direction.x.size(); i++) {
      // TODO
      //direction.active_set.bounds.at_lower_bound.erase(i);
      //direction.active_set.bounds.at_upper_bound.erase(i);
   }
   // constraints: only when p-n = 0
   for (size_t j = 0; j < direction.multipliers.constraints.size(); j++) {
      // compute constraint violation
      double constraint_violation = 0.;
      try {
         size_t i = this->elastic_variables.positive.at(j);
         constraint_violation += direction.x[i];
      }
      catch (const std::out_of_range& e) {
      }
      try {
         size_t i = this->elastic_variables.negative.at(j);
         constraint_violation += direction.x[i];
      }
      catch (const std::out_of_range& e) {
      }
      // update active set
      if (0. < constraint_violation) {
         // TODO
         //direction.active_set.constraints.at_lower_bound.erase(j);
         //direction.active_set.constraints.at_upper_bound.erase(j);
      }
   }
}