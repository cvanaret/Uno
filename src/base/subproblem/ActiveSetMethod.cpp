#include <cmath>
#include <map>
#include "ActiveSetMethod.hpp"
#include "Constraint.hpp"
#include "Vector.hpp"
#include "Logger.hpp"

ActiveSetMethod::ActiveSetMethod(Problem& problem, bool scale_residuals) :
   Subproblem(problem, L1_NORM, scale_residuals) {
}

Iterate ActiveSetMethod::evaluate_initial_point(const Problem& problem, const std::vector<double>& x, const Multipliers& multipliers) {
   Iterate first_iterate(x, multipliers);
   /* compute the optimality and feasibility measures of the initial point */
   this->compute_optimality_measures(problem, first_iterate);
   return first_iterate;
}

void ActiveSetMethod::compute_optimality_measures(const Problem& problem, Iterate& iterate) {
   // feasibility
   this->compute_residuals(problem, iterate, iterate.multipliers, 1.);
   // optimality
   iterate.compute_objective(problem);
   iterate.progress = {iterate.residuals.constraints, iterate.objective};
}

void ActiveSetMethod::compute_infeasibility_measures(const Problem& problem, Iterate& iterate, const Direction& direction) {
   iterate.compute_constraints(problem);
   // feasibility measure: residual of all constraints
   double feasibility = problem.compute_constraint_residual(iterate.constraints, this->residual_norm);
   // optimality measure: residual of linearly infeasible constraints
   double objective = problem.compute_constraint_residual(iterate.constraints, direction.constraint_partition.infeasible, this->residual_norm);
   iterate.progress = {feasibility, objective};
}

/* QP */

void ActiveSetMethod::recover_l1qp_active_set_(Problem& problem, Direction& direction, const ElasticVariables& elastic_variables) {
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
         int i = elastic_variables.positive.at(j);
         constraint_violation += direction.x[i];
      }
      catch (const std::out_of_range& e) {
      }
      try {
         int i = elastic_variables.negative.at(j);
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

void ActiveSetMethod::generate_elastic_variables_(Problem& problem, ElasticVariables& elastic_variables) {
   // generate elastic variables p and n on the fly to relax the constraints
   int elastic_index = problem.number_variables;
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (-INFINITY < problem.constraint_bounds[j].lb) {
         // nonpositive variable n that captures the negative part of the constraint violation
         elastic_variables.negative[j] = elastic_index;
         elastic_index++;
      }
      if (problem.constraint_bounds[j].ub < INFINITY) {
         // nonnegative variable p that captures the positive part of the constraint violation
         elastic_variables.positive[j] = elastic_index;
         elastic_index++;
      }
   }
}

void ActiveSetMethod::compute_l1_linear_objective_(Iterate& current_iterate, ConstraintPartition& constraint_partition) {
   /* objective function: sum of gradients of infeasible constraints */
   SparseVector objective_gradient;
   for (int j: constraint_partition.infeasible) {
      for (const auto[i, derivative]: current_iterate.constraints_jacobian[j]) {
         if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
            objective_gradient[i] -= derivative;
         }
         else {
            objective_gradient[i] += derivative;
         }
      }
   }
   current_iterate.set_objective_gradient(objective_gradient);
}

void ActiveSetMethod::generate_l1_multipliers_(Problem& problem, ConstraintPartition& constraint_partition) {
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
         this->constraints_multipliers[j] = 1.;
      }
      else if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
         this->constraints_multipliers[j] = -1.;
      }
      // otherwise, leave the multiplier as it is
   }
}

void ActiveSetMethod::generate_feasibility_bounds_(Problem& problem, std::vector<double>& current_constraints, ConstraintPartition& constraint_partition) {
   for (size_t j = 0; j < problem.number_constraints; j++) {
      double lb, ub;
      if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
         lb = -INFINITY;
         ub = problem.constraint_bounds[j].lb - current_constraints[j];
      }
      else if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
         lb = problem.constraint_bounds[j].ub - current_constraints[j];
         ub = INFINITY;
      }
      else { // FEASIBLE
         lb = problem.constraint_bounds[j].lb - current_constraints[j];
         ub = problem.constraint_bounds[j].ub - current_constraints[j];
      }
      this->constraints_bounds[j] = {lb, ub};
   }
}

/* LP */

Direction ActiveSetMethod::compute_lp_step_(Problem& problem, QPSolver& solver, Iterate& current_iterate, double trust_region_radius) {
   DEBUG << "Current point: ";
   print_vector(DEBUG, current_iterate.x);
   DEBUG << "Current constraint multipliers: ";
   print_vector(DEBUG, current_iterate.multipliers.constraints);
   DEBUG << "Current lb multipliers: ";
   print_vector(DEBUG, current_iterate.multipliers.lower_bounds);
   DEBUG << "Current ub multipliers: ";
   print_vector(DEBUG, current_iterate.multipliers.upper_bounds);

   /* bounds of the variables */
   this->generate_variables_bounds_(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->generate_constraints_bounds(problem, current_iterate.constraints);

   /* generate the initial point */
   clear(this->initial_point);

   /* solve the QP */
   Direction direction =
         solver.solve_LP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian,
               this->initial_point);
   direction.objective_multiplier = problem.objective_sign;
   direction.phase = OPTIMALITY;
   direction.predicted_reduction = [&](double step_length) {
      return ActiveSetMethod::compute_lp_predicted_reduction_(direction, step_length);
   };
   this->number_subproblems_solved++;
   DEBUG << direction;
   return direction;
}

double ActiveSetMethod::compute_lp_predicted_reduction_(Direction& direction, double step_length) {
   // the predicted reduction is linear in the step length
   return -step_length * direction.objective;
}

Direction ActiveSetMethod::compute_l1lp_step_(Problem& problem, QPSolver& solver, Iterate& current_iterate, Direction& phase_2_direction,
      double trust_region_radius) {
   DEBUG << "\nCreating the restoration problem with " << phase_2_direction.constraint_partition.infeasible.size()
         << " infeasible constraints\n";

   /* compute the objective */
   this->compute_l1_linear_objective_(current_iterate, phase_2_direction.constraint_partition);

   /* bounds of the variables */
   this->generate_variables_bounds_(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->generate_feasibility_bounds_(problem, current_iterate.constraints, phase_2_direction.constraint_partition);

   /* solve the QP */
   Direction direction =
         solver.solve_LP(this->variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate
         .constraints_jacobian,phase_2_direction.x);
   direction.objective_multiplier = 0.;
   direction.phase = RESTORATION;
   direction.constraint_partition = phase_2_direction.constraint_partition;
//    direction.predicted_reduction = [&](double step_length) {
//        return this->compute_lp_predicted_reduction_(direction, step_length);
//    };
   direction.predicted_reduction = [&](double step_length) {
      return this->compute_lp_predicted_reduction_(direction, step_length);
   };
   this->number_subproblems_solved++;
   DEBUG << direction;
   return direction;
}
