#include <cassert>
#include "Subproblem.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "optimization/Constraint.hpp"

Subproblem::Subproblem(size_t max_number_variables, size_t number_constraints, SecondOrderCorrection soc_strategy,
         bool is_second_order_method, Norm residual_norm):
      soc_strategy(soc_strategy), variable_bounds(max_number_variables),
      direction(max_number_variables, number_constraints),
      is_second_order_method(is_second_order_method), residual_norm(residual_norm) {
}

void Subproblem::initialize(Statistics& /*statistics*/, const NonlinearProblem& /*problem*/, Iterate& /*first_iterate*/) {
   // by default, do nothing
}

void Subproblem::set_variable_bounds(const NonlinearProblem& problem, const Iterate& current_iterate, double trust_region_radius) {
   // bounds intersected with trust region
   // very important: apply the trust region only on the original variables
   for (size_t i = 0; i < problem.number_variables; i++) {
      double lb = std::max(current_iterate.x[i] - trust_region_radius, problem.get_variable_lower_bound(i));
      double ub = std::min(current_iterate.x[i] + trust_region_radius, problem.get_variable_upper_bound(i));
      this->variable_bounds[i] = {lb, ub};
   }
   for (size_t i = problem.get_number_original_variables(); i < problem.number_variables; i++) {
      const double lb = problem.get_variable_lower_bound(i);
      const double ub = problem.get_variable_upper_bound(i);
      this->variable_bounds[i] = {lb, ub};
   }
}

double Subproblem::compute_first_order_error(const Model& model, Iterate& iterate) const {
   iterate.evaluate_lagrangian_gradient(model, iterate.multipliers.constraints, iterate.multipliers.lower_bounds, iterate.multipliers.upper_bounds);
   return norm(iterate.lagrangian_gradient, this->residual_norm);
}

// complementary slackness error
double Subproblem::compute_complementarity_error(const Model& model, const Iterate& iterate, const std::vector<double>& constraint_multipliers,
      const std::vector<double>& lower_bounds_multipliers, const std::vector<double>& upper_bounds_multipliers) {
   double error = 0.;
   // bound constraints
   for (size_t i = 0; i < model.number_variables; i++) {
      if (is_finite(model.get_variable_lower_bound(i))) {
         error += std::abs(lower_bounds_multipliers[i] * (iterate.x[i] - model.get_variable_lower_bound(i)));
      }
      if (is_finite(model.get_variable_upper_bound(i))) {
         error += std::abs(upper_bounds_multipliers[i] * (iterate.x[i] - model.get_variable_upper_bound(i)));
      }
   }
   // constraints
   for (size_t j = 0; j < model.number_constraints; j++) {
      const double multiplier_j = constraint_multipliers[j];
      const double lower_bound = model.get_constraint_lower_bound(j);
      const double upper_bound = model.get_constraint_upper_bound(j);
      if (iterate.original_evaluations.constraints[j] < lower_bound) { // violated lower
         // the optimal multiplier is 1
         error += std::abs((1. - multiplier_j) * (lower_bound - iterate.original_evaluations.constraints[j]));
      }
      else if (upper_bound < iterate.original_evaluations.constraints[j]) { // violated upper
         // the optimal multiplier is -1
         error += std::abs((1. + multiplier_j) * (iterate.original_evaluations.constraints[j] - upper_bound));
      }
      else if (is_finite(lower_bound) && 0. < multiplier_j) { // lower bound
         error += std::abs(multiplier_j * (iterate.original_evaluations.constraints[j] - lower_bound));
      }
      else if (is_finite(upper_bound) && multiplier_j < 0.) { // upper bound
         error += std::abs(multiplier_j * (iterate.original_evaluations.constraints[j] - upper_bound));
      }
   }
   return error;
}

void Subproblem::compute_nonlinear_residuals(const NonlinearProblem& problem, Iterate& iterate) const {
   iterate.evaluate_constraints(problem.model);
   iterate.nonlinear_errors.constraints = problem.compute_constraint_violation(iterate.original_evaluations.constraints, L1_NORM);
   iterate.nonlinear_errors.stationarity = this->compute_first_order_error(problem.model, iterate);
   iterate.nonlinear_errors.complementarity = Subproblem::compute_complementarity_error(problem.model, iterate, iterate.multipliers.constraints,
         iterate.multipliers.lower_bounds, iterate.multipliers.upper_bounds);
}

Direction Subproblem::compute_second_order_correction(const NonlinearProblem& /*problem*/, Iterate& /*trial_iterate*/) {
   assert(false && "Subproblem::compute_second_order_correction");
}

void Subproblem::postprocess_accepted_iterate(const NonlinearProblem& /*problem*/, Iterate& /*iterate*/) {
   // by default, do nothing
}