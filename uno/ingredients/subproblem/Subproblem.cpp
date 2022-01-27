#include <cassert>
#include "Subproblem.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "optimization/Constraint.hpp"

Subproblem::Subproblem(size_t number_variables, size_t max_number_variables, size_t number_constraints, SecondOrderCorrection soc_strategy,
         bool is_second_order_method, Norm residual_norm):
      number_variables(number_variables), max_number_variables(max_number_variables), number_constraints(number_constraints),
      soc_strategy(soc_strategy), current_variable_bounds(max_number_variables),
      objective_gradient(max_number_variables), // SparseVector
      constraint_bounds(number_constraints), direction(max_number_variables, number_constraints),
      is_second_order_method(is_second_order_method), residual_norm(residual_norm) {
}

void Subproblem::initialize(Statistics& /*statistics*/, const Problem& /*problem*/, Iterate& /*first_iterate*/) {
   // by default, do nothing
}

double Subproblem::get_proximal_coefficient() const {
   // by default, return 1.
   return 1.;
}

double Subproblem::push_variable_to_interior(double variable_value, const Range& variable_bounds) {
   const double k1 = 1e-2;
   const double k2 = 1e-2;

   const double range = variable_bounds.ub - variable_bounds.lb;
   const double perturbation_lb = std::min(k1 * std::max(1., std::abs(variable_bounds.lb)), k2 * range);
   const double perturbation_ub = std::min(k1 * std::max(1., std::abs(variable_bounds.ub)), k2 * range);
   variable_value = std::max(variable_value, variable_bounds.lb + perturbation_lb);
   variable_value = std::min(variable_value, variable_bounds.ub - perturbation_ub);
   return variable_value;
}

void Subproblem::set_current_variable_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) {
   // bounds intersected with trust region
   // very important: apply the trust region only on the original variables
   for (size_t i = 0; i < problem.number_variables; i++) {
      const double lb = std::max(-trust_region_radius, problem.get_variable_lower_bound(i) - current_iterate.x[i]);
      const double ub = std::min(trust_region_radius, problem.get_variable_upper_bound(i) - current_iterate.x[i]);
      this->current_variable_bounds[i] = {lb, ub};
   }
}

void Subproblem::set_constraint_bounds(const Problem& problem, const std::vector<double>& current_constraints) {
   for (size_t j = 0; j < problem.number_constraints; j++) {
      const double lb = problem.get_constraint_lower_bound(j) - current_constraints[j];
      const double ub = problem.get_constraint_upper_bound(j) - current_constraints[j];
      this->constraint_bounds[j] = {lb, ub};
   }
}

void Subproblem::set_scaled_objective_gradient(const Problem& problem, Iterate& current_iterate, double objective_multiplier) {
   // scale objective gradient
   if (objective_multiplier == 0.) {
      this->objective_gradient.clear();
   }
   else {
      current_iterate.evaluate_objective_gradient(problem);
      this->objective_gradient = current_iterate.objective_gradient;
      if (objective_multiplier != 1.) {
         scale(this->objective_gradient, objective_multiplier);
      }
   }
}

double Subproblem::compute_first_order_error(const Problem& problem, Iterate& iterate, double objective_multiplier) const {
   iterate.evaluate_lagrangian_gradient(problem, objective_multiplier, iterate.multipliers.constraints, iterate.multipliers.lower_bounds,
         iterate.multipliers.upper_bounds);
   return norm(iterate.lagrangian_gradient, this->residual_norm);
}

// complementary slackness error
double Subproblem::compute_complementarity_error(const Problem& problem, const Iterate& iterate, const std::vector<double>& constraint_multipliers,
      const std::vector<double>& lower_bounds_multipliers, const std::vector<double>& upper_bounds_multipliers) {
   double error = 0.;
   // bound constraints
   for (size_t i = 0; i < problem.number_variables; i++) {
      if (is_finite_lower_bound(problem.get_variable_lower_bound(i))) {
         error += std::abs(lower_bounds_multipliers[i] * (iterate.x[i] - problem.get_variable_lower_bound(i)));
      }
      if (is_finite_upper_bound(problem.get_variable_upper_bound(i))) {
         error += std::abs(upper_bounds_multipliers[i] * (iterate.x[i] - problem.get_variable_upper_bound(i)));
      }
   }
   // constraints
   for (size_t j = 0; j < problem.number_constraints; j++) {
      const double multiplier_j = constraint_multipliers[j];
      const double lower_bound = problem.get_constraint_lower_bound(j);
      const double upper_bound = problem.get_constraint_upper_bound(j);
      if (iterate.constraints[j] < lower_bound) {
         // violated lower: the multiplier is 1 at optimum
         error += std::abs((1. - multiplier_j) * (lower_bound - iterate.constraints[j]));
      }
      else if (upper_bound < iterate.constraints[j]) {
         // violated upper: the multiplier is -1 at optimum
         error += std::abs((1. + multiplier_j) * (iterate.constraints[j] - upper_bound));
      }
      else if (is_finite_lower_bound(lower_bound) && 0. < multiplier_j) {
         error += std::abs(multiplier_j * (iterate.constraints[j] - lower_bound));
      }
      else if (is_finite_upper_bound(upper_bound) && multiplier_j < 0.) {
         error += std::abs(multiplier_j * (iterate.constraints[j] - upper_bound));
      }
   }
   return error;
}

void Subproblem::compute_progress_measures(const Problem& problem, Iterate& iterate) {
   // feasibility measure: constraint violation
   iterate.evaluate_constraints(problem);
   iterate.nonlinear_errors.constraints = problem.compute_constraint_violation(iterate.x, iterate.constraints);
   // optimality measure: objective value
   iterate.evaluate_objective(problem);
   iterate.progress = {iterate.nonlinear_errors.constraints, iterate.objective};
}

void Subproblem::compute_residuals(const Problem& problem, Iterate& iterate, double objective_multiplier) const {
   iterate.evaluate_constraints(problem);
   iterate.nonlinear_errors.constraints = problem.compute_constraint_violation(iterate.x, iterate.constraints);
   // compute the KKT error only if the objective multiplier is positive
   iterate.nonlinear_errors.KKT = this->compute_first_order_error(problem, iterate, 0. < objective_multiplier ? objective_multiplier : 1.);
   iterate.nonlinear_errors.FJ = this->compute_first_order_error(problem, iterate, 0.);
   iterate.nonlinear_errors.complementarity = Subproblem::compute_complementarity_error(problem, iterate, iterate.multipliers.constraints,
         iterate.multipliers.lower_bounds, iterate.multipliers.upper_bounds);
}

Direction Subproblem::compute_second_order_correction(const Problem& /*problem*/, Iterate& /*trial_iterate*/) {
   assert(false && "Subproblem::compute_second_order_correction");
}

void Subproblem::register_accepted_iterate(Iterate& /*iterate*/) {
   // by default, do nothing
}