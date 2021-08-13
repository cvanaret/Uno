#include <cassert>
#include "Subproblem.hpp"

Subproblem::Subproblem(size_t number_variables, size_t number_constraints) :
      number_variables(number_variables),
      number_constraints(number_constraints),
      variables_bounds(number_variables),
      constraints_multipliers(number_constraints),
      // objective_gradient is a SparseVector
      // constraints_jacobian is a vector of SparseVectors
      constraints_jacobian(number_constraints),
      constraints_bounds(number_constraints),
      number_subproblems_solved(0), subproblem_definition_changed(false) {
}

void Subproblem::evaluate_constraints(const Problem& problem, Iterate& iterate) const {
   iterate.compute_constraints(problem);
}

Iterate Subproblem::generate_initial_iterate(Statistics& /*statistics*/, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
   Iterate first_iterate(x, multipliers);
   /* compute the optimality and feasibility measures of the initial point */
   this->evaluate_constraints(problem, first_iterate);
   this->compute_progress_measures(problem, first_iterate);
   return first_iterate;
}

void Subproblem::compute_progress_measures(const Problem& problem, Iterate& iterate) {
   iterate.compute_constraints(problem);
   // feasibility measure: residual of all constraints
   iterate.errors.constraints = problem.compute_constraint_violation(iterate.constraints, L1_NORM);
   // optimality
   iterate.compute_objective(problem);
   iterate.progress = {iterate.errors.constraints, iterate.objective};
}

double Subproblem::push_variable_to_interior(double variable_value, const Range& variable_bounds) {
   double k1 = 1e-2;
   double k2 = 1e-2;

   double perturbation_lb = std::min(k1 * std::max(1., std::abs(variable_bounds.lb)), k2 * (variable_bounds.ub - variable_bounds.lb));
   double perturbation_ub = std::min(k1 * std::max(1., std::abs(variable_bounds.ub)), k2 * (variable_bounds.ub - variable_bounds.lb));
   variable_value = std::max(variable_value, variable_bounds.lb + perturbation_lb);
   variable_value = std::min(variable_value, variable_bounds.ub - perturbation_ub);
   return variable_value;
}

void Subproblem::set_trust_region(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) {
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);
}

void Subproblem::set_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) {
   /* bounds intersected with trust region  */
   // very important: apply the trust region only on the original variables
   for (size_t i = 0; i < problem.number_variables; i++) {
      double lb = std::max(-trust_region_radius, problem.variables_bounds[i].lb - current_iterate.x[i]);
      double ub = std::min(trust_region_radius, problem.variables_bounds[i].ub - current_iterate.x[i]);
      this->variables_bounds[i] = {lb, ub};
   }
}

void Subproblem::set_constraints_bounds(const Problem& problem, const std::vector<double>& current_constraints) {
   for (size_t j = 0; j < problem.number_constraints; j++) {
      double lb = problem.constraint_bounds[j].lb - current_constraints[j];
      double ub = problem.constraint_bounds[j].ub - current_constraints[j];
      this->constraints_bounds[j] = {lb, ub};
   }
}

void Subproblem::compute_feasibility_linear_objective(const Iterate& current_iterate, const ConstraintPartition& constraint_partition) {
   /* objective function: sum of gradients of infeasible constraints */
   this->objective_gradient.clear();
   for (int j: constraint_partition.infeasible) {
      for (const auto [i, derivative]: current_iterate.constraints_jacobian[j]) {
         if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
            this->objective_gradient[i] -= derivative;
         }
         else {
            this->objective_gradient[i] += derivative;
         }
      }
   }
}

void Subproblem::generate_feasibility_bounds(const Problem& problem, const std::vector<double>& current_constraints, const ConstraintPartition&
constraint_partition) {
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

double Subproblem::compute_first_order_error(const Problem& problem, Iterate& iterate, double objective_multiplier) {
   std::vector<double> lagrangian_gradient = iterate.lagrangian_gradient(problem, objective_multiplier, iterate.multipliers);
   return norm_1(lagrangian_gradient);
}

/* complementary slackness error. Use abs/1e-8 to safeguard */
double Subproblem::compute_complementarity_error(const Problem& problem, Iterate& iterate, const Multipliers& multipliers) {
   double error = 0.;
   /* bound constraints */
   for (size_t i = 0; i < problem.number_variables; i++) {
      if (-INFINITY < problem.variables_bounds[i].lb) {
         error += std::abs(multipliers.lower_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].lb));
      }
      if (problem.variables_bounds[i].ub < INFINITY) {
         error += std::abs(multipliers.upper_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].ub));
      }
   }
   /* constraints */
   iterate.compute_constraints(problem);
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

double Subproblem::compute_constraint_violation(const Problem& problem, const Iterate& iterate) const {
   return problem.compute_constraint_violation(iterate.constraints, L1_NORM);
}

void Subproblem::compute_errors(const Problem& problem, Iterate& iterate, double objective_multiplier) const {
   iterate.compute_constraints(problem);
   iterate.errors.constraints = this->compute_constraint_violation(problem, iterate);
   // compute the KKT error only if the objective multiplier is positive
   iterate.errors.KKT = Subproblem::compute_first_order_error(problem, iterate, 0. < objective_multiplier ? objective_multiplier : 1.);
   iterate.errors.FJ = Subproblem::compute_first_order_error(problem, iterate, 0.);
   iterate.errors.complementarity = Subproblem::compute_complementarity_error(problem, iterate, iterate.multipliers);
}

Direction Subproblem::compute_second_order_correction(const Problem& /*problem*/, Iterate& /*trial_iterate*/) {
   assert(false && "Subproblem::compute_second_order_correction");
}

void Subproblem::register_accepted_iterate(Iterate& /*iterate*/) {
   // by default, do nothing
}