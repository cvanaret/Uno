#include <cmath>
#include <map>
#include <memory>
#include "Sl1QP.hpp"
#include "Constraint.hpp"
#include "Vector.hpp"
#include "BQPDSolver.hpp"
#include "QPSolverFactory.hpp"

Sl1QP::Sl1QP(const Problem& problem, const std::string& QP_solver, const std::string& hessian_evaluation_method, bool use_trust_region, bool
scale_residuals, double initial_parameter) :
// compute the number of variables and call the private constructor
      Sl1QP(problem, QP_solver, hessian_evaluation_method, use_trust_region, scale_residuals, initial_parameter,
            this->count_elastic_variables_(problem)) {
}

size_t Sl1QP::count_elastic_variables_(const Problem& problem) {
   size_t number_variables = problem.number_variables;
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

Sl1QP::Sl1QP(const Problem& problem, const std::string& QP_solver, const std::string& hessian_evaluation_method, bool use_trust_region, bool
scale_residuals, double initial_parameter, int number_variables) : ActiveSetMethod(problem, scale_residuals), solver(
      QPSolverFactory::create(QP_solver, number_variables, problem.number_constraints,
            problem.hessian_maximum_number_nonzeros + problem.number_variables, true)),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      hessian_evaluation(HessianEvaluationFactory::create(hessian_evaluation_method,
            problem.number_variables, problem.hessian_maximum_number_nonzeros, !use_trust_region)),
      penalty_parameter(initial_parameter), parameters({10., 0.1, 0.1}),
      number_variables(number_variables) {
   // generate elastic variables p and n on the fly to relax the constraints
   this->generate_elastic_variables_(problem, this->elastic_variables_);
   // TODO let the solver resize the Hessian space
}

void Sl1QP::generate(const Problem& /*problem*/, const Iterate& /*current_iterate*/, double /*objective_multiplier*/, double
/*trust_region_radius*/) {
}

void Sl1QP::update_objective_multipliers(const Problem& /*problem*/, const Iterate& /*current_iterate*/, double /*objective_multiplier*/) {
}

Direction Sl1QP::compute_direction(const Problem& problem, Iterate& current_iterate, double trust_region_radius) {
   DEBUG << "penalty parameter: " << this->penalty_parameter << "\n";

   // evaluate constraints
   current_iterate.compute_constraints(problem);

   /* stage a: compute the step within trust region */
   Direction direction = this->solve_l1qp_subproblem_(problem, current_iterate, trust_region_radius, this->penalty_parameter);

   /* penalty update: if penalty parameter is already 0, no need to decrease it */
   if (0. < this->penalty_parameter) {
      /* check infeasibility */
      double linearized_residual = this->compute_linearized_constraint_residual_(direction.x);
      DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";

      // if problem had to be relaxed
      if (linearized_residual != 0.) {
         double current_penalty_parameter = this->penalty_parameter;

         /* stage c: solve the ideal l1 penalty problem with a zero penalty (no objective) */
         DEBUG << "Compute ideal solution:\n";
         Direction ideal_direction = this->solve_l1qp_subproblem_(problem, current_iterate, trust_region_radius, 0.);

         /* compute the ideal error (with a zero penalty parameter) */
         double ideal_error = this->compute_error_(problem, current_iterate, ideal_direction.multipliers, 0.);
         DEBUG << "Ideal error: " << ideal_error << "\n";

         if (ideal_error == 0.) {
            /* stage f: update the penalty parameter */
            this->penalty_parameter = 0.;
         }
         else {
            double ideal_linearized_residual = this->compute_linearized_constraint_residual_(ideal_direction.x);
            DEBUG << "Linearized residual mk(dk): " << ideal_linearized_residual << "\n\n";

            /* decrease penalty parameter to satisfy 2 conditions */
            bool condition1 = false, condition2 = false;
            while (!condition2) {
               if (!condition1) {
                  /* stage d: reach a fraction of the ideal decrease */
                  if ((ideal_linearized_residual == 0. && linearized_residual == 0) || (ideal_linearized_residual != 0. &&
                                                                                        current_iterate.progress.feasibility -
                                                                                        linearized_residual >= this->parameters.epsilon1 *
                                                                                                               (current_iterate
                                                                                                                      .progress.feasibility -
                                                                                                                ideal_linearized_residual))) {
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
                     direction = this->solve_l1qp_subproblem_(problem, current_iterate, trust_region_radius, this->penalty_parameter);

                     linearized_residual = this->compute_linearized_constraint_residual_(direction.x);
                     DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";
                  }
               }
            }

            /* stage f: update the penalty parameter */
            double updated_penalty_parameter = this->penalty_parameter;
            double term = ideal_error / std::max(1., current_iterate.progress.feasibility);
            this->penalty_parameter = std::min(this->penalty_parameter, term * term);
            if (this->penalty_parameter < updated_penalty_parameter) {
               direction = this->solve_l1qp_subproblem_(problem, current_iterate, trust_region_radius, this->penalty_parameter);
            }

         }

         if (this->penalty_parameter < current_penalty_parameter) {
            DEBUG << "\n*** Penalty parameter updated to " << this->penalty_parameter << "\n";
            this->subproblem_definition_changed = true;
            if (this->penalty_parameter == 0.) {
               direction = ideal_direction;
            }
         }
      }
   }

   direction.objective_multiplier = penalty_parameter;
   direction.predicted_reduction = [&](double step_length) {
      return this->compute_predicted_reduction_(problem, current_iterate, direction, step_length);
   };
   return direction;
}

Direction Sl1QP::compute_l1qp_step_(const Problem& problem, QPSolver& solver, Iterate& current_iterate, ConstraintPartition& constraint_partition,
      std::vector<double>& initial_solution, double trust_region_radius) {
   /* compute the objective */
   this->compute_l1_linear_objective_(current_iterate, constraint_partition);

   /* bounds of the variables */
   this->generate_variables_bounds_(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->generate_feasibility_bounds_(problem, current_iterate.constraints, constraint_partition);

   /* solve the QP */
   Direction direction =
         solver.solve_QP(this->variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate
         .constraints_jacobian,this->hessian_evaluation->hessian, initial_solution);
   direction.objective_multiplier = 0.;
   direction.constraint_partition = constraint_partition;
   this->number_subproblems_solved++;
   DEBUG << direction;
   return direction;
}

Direction Sl1QP::compute_l1qp_step_(const Problem& problem, QPSolver& solver, Iterate& current_iterate, double penalty_parameter,
      ElasticVariables& elastic_variables, double trust_region_radius) {
   current_iterate.compute_objective_gradient(problem);
   SparseVector objective_gradient;
   if (penalty_parameter != 0.) {
      for (const auto[i, derivative]: current_iterate.objective_gradient) {
         objective_gradient[i] = penalty_parameter * derivative;
      }
   }
   /* add contribution of positive part variables */
   for (const auto[j, i]: elastic_variables.positive) {
      current_iterate.constraints_jacobian[j][i] = -1.;
      objective_gradient[i] = 1.;
   }
   /* add contribution of negative part variables */
   for (const auto[j, i]: elastic_variables.negative) {
      current_iterate.constraints_jacobian[j][i] = 1.;
      objective_gradient[i] = 1.;
   }
   //current_iterate.set_objective_gradient(objective_gradient);

   /* bounds of the variables */
   this->generate_variables_bounds_(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->generate_constraints_bounds(problem, current_iterate.constraints);

   /* generate the initial point */
   for (size_t i = 0; i < problem.number_variables; i++) {
      this->initial_point[i] = 0.;
   }

   DEBUG << "Bounds:\n";
   for (size_t i = 0; i < this->variables_bounds.size(); i++) {
      DEBUG << "x" << i << " in [" << this->variables_bounds[i].lb << ", " << this->variables_bounds[i].ub << "]\n";
   }
   DEBUG << "Hessian: " << this->hessian_evaluation->hessian << "\n";
   DEBUG << "Objective gradient: ";
   print_vector(DEBUG, current_iterate.objective_gradient);
   for (size_t j = 0; j < constraints_bounds.size(); j++) {
      DEBUG << "Constraint " << j << ": ";
      print_vector(DEBUG, current_iterate.constraints_jacobian[j], ' ');
      DEBUG << " in [" << constraints_bounds[j].lb << ", " << constraints_bounds[j].ub << "]\n";
   }

   /* solve the QP */
   Direction direction = solver.solve_QP(this->variables_bounds, this->constraints_bounds, objective_gradient,
         current_iterate.constraints_jacobian, this->hessian_evaluation->hessian, this->initial_point);
   direction.phase = OPTIMALITY;
   this->number_subproblems_solved++;
   // recompute active set: constraints are active when p-n = 0
   this->recover_l1qp_active_set_(problem, direction, elastic_variables);
   DEBUG << direction;

   /* remove p and n */
   direction.x.resize(current_iterate.x.size());
   direction.multipliers.lower_bounds.resize(current_iterate.x.size());
   direction.multipliers.upper_bounds.resize(current_iterate.x.size());
   direction.norm = norm_inf(direction.x);

   /* remove contribution of positive part variables */
   for (const auto[j, i]: elastic_variables.positive) {
      current_iterate.constraints_jacobian[j].erase(i);
      current_iterate.objective_gradient.erase(i);
   }
   /* remove contribution of negative part variables */
   for (const auto[j, i]: elastic_variables.negative) {
      current_iterate.constraints_jacobian[j].erase(i);
      current_iterate.objective_gradient.erase(i);
   }
   return direction;
}

Direction Sl1QP::solve_l1qp_subproblem_(const Problem& problem, Iterate& current_iterate, double trust_region_radius, double penalty_parameter) {
   /* compute l1QP step */
   this->hessian_evaluation->compute(problem, current_iterate.x, penalty_parameter, current_iterate.multipliers.constraints);
   Direction direction = this->compute_l1qp_step_(problem, *this->solver, current_iterate, penalty_parameter, this->elastic_variables_,
         trust_region_radius);
   return direction;
}

Direction Sl1QP::restore_feasibility(const Problem&, Iterate&, Direction&, double) {
   throw std::out_of_range("Sl1QP.compute_infeasibility_step is not implemented, since l1QP are always feasible");
}

double Sl1QP::compute_predicted_reduction_(const Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) {
   // the predicted reduction is quadratic
   if (step_length == 1.) {
      return current_iterate.progress.feasibility - direction.objective;
   }
   else {
      double linear_term = dot(direction.x, current_iterate.objective_gradient);
      double quadratic_term = this->hessian_evaluation->hessian.quadratic_product(direction.x, direction.x) / 2.;
      // determine the constraint violation term: c(x_k) + alpha*\nabla c(x_k)^T d
      std::vector<double> scaled_constraints(current_iterate.constraints);
      for (size_t j = 0; j < current_iterate.constraints.size(); j++) {
         scaled_constraints[j] += step_length * dot(direction.x, current_iterate.constraints_jacobian[j]);
      }
      double constraint_violation = problem.compute_constraint_residual(scaled_constraints, this->residual_norm);
      return current_iterate.progress.feasibility - constraint_violation - step_length * (linear_term + step_length * quadratic_term);
   }
}

/* private methods */

void Sl1QP::generate_variables_bounds_(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) {
   // p and n are nonnegative
   //this->variables_bounds(this->number_variables, {0., INFINITY});

   /* original bounds intersected with trust region  */
   for (size_t i = 0; i < problem.number_variables; i++) {
      double lb = std::max(-trust_region_radius, this->variables_bounds[i].lb - current_iterate.x[i]);
      double ub = std::min(trust_region_radius, this->variables_bounds[i].ub - current_iterate.x[i]);
      this->variables_bounds[i] = {lb, ub};
   }
}

double Sl1QP::compute_linearized_constraint_residual_(std::vector<double>& direction) {
   double residual = 0.;
   // l1 residual of the linearized constraints
   for (std::pair<const int, int>& element: this->elastic_variables_.positive) {
      int i = element.second;
      residual += direction[i];
   }
   for (std::pair<const int, int>& element: this->elastic_variables_.negative) {
      int i = element.second;
      residual += direction[i];
   }
   return residual;
}

double Sl1QP::compute_error_(const Problem& problem, Iterate& iterate, Multipliers& multipliers, double penalty_parameter) {
   /* measure that combines KKT error and complementarity error */
   double error = 0.;

   /* KKT error */
   std::vector<double> lagrangian_gradient = iterate.lagrangian_gradient(problem, penalty_parameter, multipliers);
   error += norm_1(lagrangian_gradient);
   /* complementarity error */
   error += this->compute_complementarity_error_(problem, iterate, multipliers);
   return error;
}

/* complementary slackness error. Use abs/1e-8 to safeguard */
double Sl1QP::compute_complementarity_error_(const Problem& problem, Iterate& iterate, const Multipliers& multipliers) const {
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
