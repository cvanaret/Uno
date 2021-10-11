#include <cassert>
#include "Subproblem.hpp"

Subproblem::Subproblem(size_t number_variables, size_t max_number_variables, size_t number_constraints, SecondOrderCorrection soc_strategy) :
      number_variables(number_variables), max_number_variables(max_number_variables), number_constraints(number_constraints),
      soc_strategy(soc_strategy), variables_bounds(max_number_variables), constraints_multipliers(number_constraints),
      objective_gradient(max_number_variables), // SparseVector
      constraint_jacobian(number_constraints), // vector of SparseVectors
      constraints_bounds(number_constraints), direction(max_number_variables, number_constraints) {
   for (auto& constraint_gradient: this->constraint_jacobian) {
      constraint_gradient.reserve(this->max_number_variables);
   }
}

void Subproblem::initialize(Statistics& /*statistics*/, const Problem& problem, Iterate& first_iterate) {
   /* compute the optimality and feasibility measures of the initial point */
   first_iterate.evaluate_constraints(problem);
   this->compute_progress_measures(problem, first_iterate);
}

void Subproblem::add_variable(size_t i, double /*current_value*/, const Range& bounds, double objective_term, size_t j, double jacobian_term) {
   assert(i < this->max_number_variables && "The index is larger than the preallocated size");
   assert(j < this->number_constraints && "The constraint index is larger than the preallocated size");
   this->variables_bounds[i] = bounds;
   this->objective_gradient.insert(i, objective_term);
   this->constraint_jacobian[j].insert(i, jacobian_term);
   this->number_variables++;
}

void Subproblem::remove_variable(size_t i, size_t j) {
   assert(i < this->max_number_variables && "The variable index is larger than the preallocated size");
   assert(j < this->number_constraints && "The constraint index is larger than the preallocated size");
   this->objective_gradient.erase(i);
   this->constraint_jacobian[j].erase(i);
   this->number_variables--;
}

void Subproblem::compute_progress_measures(const Problem& problem, Iterate& iterate) {
   iterate.evaluate_constraints(problem);
   // feasibility measure: residual of all constraints
   iterate.errors.constraints = problem.compute_constraint_violation(iterate.constraints, L1_NORM);
   // optimality
   iterate.evaluate_objective(problem);
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

void Subproblem::set_scaled_objective_gradient(const Problem& problem, Iterate& current_iterate, double objective_multiplier) {
   // scale objective gradient
   current_iterate.evaluate_objective_gradient(problem);
   if (objective_multiplier == 0.) {
      this->objective_gradient.clear();
   }
   else {
      this->objective_gradient = current_iterate.objective_gradient;
      if (objective_multiplier != 1.) {
         scale(this->objective_gradient, objective_multiplier);
      }
   }
}

void Subproblem::compute_feasibility_linear_objective(const Iterate& current_iterate, const ConstraintPartition& constraint_partition) {
   /* objective function: sum of gradients of infeasible constraints */
   this->objective_gradient.clear();
   for (size_t j: constraint_partition.lower_bound_infeasible) {
      current_iterate.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
         this->objective_gradient.insert(i, -derivative);
      });
   }
   for (size_t j: constraint_partition.upper_bound_infeasible) {
      current_iterate.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
         this->objective_gradient.insert(i, derivative);
      });
   }
}

void Subproblem::generate_feasibility_bounds(const Problem& problem, const std::vector<double>& current_constraints, const ConstraintPartition&
constraint_partition) {
   for (size_t j: constraint_partition.feasible) {
      this->constraints_bounds[j] = {problem.constraint_bounds[j].lb - current_constraints[j], problem.constraint_bounds[j].ub - current_constraints[j]};
   }
   for (size_t j: constraint_partition.lower_bound_infeasible) {
      this->constraints_bounds[j] = {-std::numeric_limits<double>::infinity(), problem.constraint_bounds[j].lb - current_constraints[j]};
   }
   for (size_t j: constraint_partition.upper_bound_infeasible) {
      this->constraints_bounds[j] = {problem.constraint_bounds[j].ub - current_constraints[j], std::numeric_limits<double>::infinity()};
   }
}

double Subproblem::compute_first_order_error(const Problem& problem, Iterate& iterate, double objective_multiplier) {
   iterate.evaluate_lagrangian_gradient(problem, objective_multiplier, iterate.multipliers);
   return norm_1(iterate.lagrangian_gradient);
}

/* complementary slackness error */
double Subproblem::compute_complementarity_error(const Problem& problem, Iterate& iterate, const Multipliers& multipliers) {
   double error = 0.;
   /* bound constraints */
   for (size_t i = 0; i < problem.number_variables; i++) {
      if (-std::numeric_limits<double>::infinity() < problem.variables_bounds[i].lb) {
         error += std::abs(multipliers.lower_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].lb));
      }
      if (problem.variables_bounds[i].ub < std::numeric_limits<double>::infinity()) {
         error += std::abs(multipliers.upper_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].ub));
      }
   }
   /* constraints */
   iterate.evaluate_constraints(problem);
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
      else if (-std::numeric_limits<double>::infinity() < problem.constraint_bounds[j].lb && 0. < multiplier_j) {
         error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraint_bounds[j].lb));
      }
      else if (problem.constraint_bounds[j].ub < std::numeric_limits<double>::infinity() && multiplier_j < 0.) {
         error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraint_bounds[j].ub));
      }
   }
   return error;
}

void Subproblem::compute_optimality_conditions(const Problem& problem, Iterate& iterate, double objective_multiplier) const {
   iterate.evaluate_objective(problem);
   iterate.evaluate_constraints(problem);
   iterate.errors.constraints = problem.compute_constraint_violation(iterate.constraints, L1_NORM);
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