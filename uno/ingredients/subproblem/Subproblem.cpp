#include <cassert>
#include "Subproblem.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "optimization/Constraint.hpp"

Subproblem::Subproblem(size_t number_variables, size_t max_number_variables, size_t number_constraints, bool uses_slacks,
      SecondOrderCorrection soc_strategy, bool is_second_order_method, Norm residual_norm) :
      number_variables(number_variables), max_number_variables(max_number_variables), number_constraints(number_constraints),
      uses_slacks(uses_slacks),
      soc_strategy(soc_strategy), variables_bounds(max_number_variables), constraints_multipliers(number_constraints),
      objective_gradient(max_number_variables), // SparseVector
      constraints_jacobian(number_constraints), // vector of SparseVectors
      constraints_bounds(number_constraints), direction(max_number_variables, number_constraints),
      is_second_order_method(is_second_order_method), residual_norm(residual_norm) {
   for (auto& constraint_gradient: this->constraints_jacobian) {
      constraint_gradient.reserve(this->max_number_variables);
   }
}

void Subproblem::initialize(Statistics& /*statistics*/, const Problem& problem, const Scaling& scaling, Iterate& first_iterate) {
   // compute the optimality and feasibility measures of the initial point
   first_iterate.evaluate_constraints(problem, scaling);
   this->compute_progress_measures(problem, scaling, first_iterate);
}

void Subproblem::add_elastic_variable(size_t i, double objective_term, size_t j, double jacobian_term) {
   assert(i < this->max_number_variables && "The variable index is larger than the preallocated size");
   assert(j < this->number_constraints && "The constraint index is larger than the preallocated size");
   this->variables_bounds[i] = {0., std::numeric_limits<double>::infinity()};
   this->objective_gradient.insert(i, objective_term);
   this->constraints_jacobian[j].insert(i, jacobian_term);
   this->number_variables++;
}

void Subproblem::remove_elastic_variable(size_t i, size_t j) {
   assert(i < this->max_number_variables && "The variable index is larger than the preallocated size");
   assert(j < this->number_constraints && "The constraint index is larger than the preallocated size");
   this->objective_gradient.erase(i);
   this->constraints_jacobian[j].erase(i);
   this->number_variables--;
}

void Subproblem::compute_progress_measures(const Problem& problem, const Scaling& scaling, Iterate& iterate) {
   iterate.evaluate_constraints(problem, scaling);
   // feasibility measure: residual of all constraints
   iterate.errors.constraints = problem.compute_constraint_violation(scaling, iterate.constraints, this->residual_norm);
   // optimality
   iterate.evaluate_objective(problem, scaling);
   iterate.progress = {iterate.errors.constraints, iterate.objective};
}

double Subproblem::push_variable_to_interior(double variable_value, const Range& variable_bounds) {
   const double k1 = 1e-2;
   const double k2 = 1e-2;

   const double perturbation_lb = std::min(k1 * std::max(1., std::abs(variable_bounds.lb)), k2 * (variable_bounds.ub - variable_bounds.lb));
   const double perturbation_ub = std::min(k1 * std::max(1., std::abs(variable_bounds.ub)), k2 * (variable_bounds.ub - variable_bounds.lb));
   variable_value = std::max(variable_value, variable_bounds.lb + perturbation_lb);
   variable_value = std::min(variable_value, variable_bounds.ub - perturbation_ub);
   return variable_value;
}

void Subproblem::set_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) {
   // bounds intersected with trust region
   // very important: apply the trust region only on the original variables
   for (size_t i = 0; i < problem.number_variables; i++) {
      const double lb = std::max(-trust_region_radius, problem.variables_bounds[i].lb - current_iterate.x[i]);
      const double ub = std::min(trust_region_radius, problem.variables_bounds[i].ub - current_iterate.x[i]);
      this->variables_bounds[i] = {lb, ub};
   }
}

void Subproblem::set_constraints_bounds(const Problem& problem, const std::vector<double>& current_constraints) {
   for (size_t j = 0; j < problem.number_constraints; j++) {
      const double lb = problem.constraint_bounds[j].lb - current_constraints[j];
      const double ub = problem.constraint_bounds[j].ub - current_constraints[j];
      this->constraints_bounds[j] = {lb, ub};
   }
}

void Subproblem::set_scaled_objective_gradient(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, double objective_multiplier) {
   // scale objective gradient
   if (objective_multiplier == 0.) {
      this->objective_gradient.clear();
   }
   else {
      current_iterate.evaluate_objective_gradient(problem, scaling);
      this->objective_gradient = current_iterate.objective_gradient;
      if (objective_multiplier != 1.) {
         scale(this->objective_gradient, objective_multiplier);
      }
   }
}

void Subproblem::compute_feasibility_linear_objective(const Iterate& current_iterate, const ConstraintPartition& constraint_partition) {
   // objective function: sum of gradients of infeasible constraints
   this->objective_gradient.clear();
   for (size_t j: constraint_partition.lower_bound_infeasible) {
      current_iterate.constraints_jacobian[j].for_each([&](size_t i, double derivative) {
         this->objective_gradient.insert(i, -derivative);
      });
   }
   for (size_t j: constraint_partition.upper_bound_infeasible) {
      current_iterate.constraints_jacobian[j].for_each([&](size_t i, double derivative) {
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

double Subproblem::compute_first_order_error(const Problem& problem, const Scaling& scaling, Iterate& iterate, double objective_multiplier) const {
   iterate.evaluate_lagrangian_gradient(problem, scaling, objective_multiplier, iterate.multipliers.constraints, iterate.multipliers.lower_bounds,
         iterate.multipliers.upper_bounds);
   return norm(iterate.lagrangian_gradient, this->residual_norm);
}

// complementary slackness error
double Subproblem::compute_complementarity_error(const Problem& problem, const Scaling& scaling, Iterate& iterate,
      const std::vector<double>& constraint_multipliers, const std::vector<double>& lower_bounds_multipliers,
      const std::vector<double>& upper_bounds_multipliers) {
   double error = 0.;
   // bound constraints
   for (size_t i = 0; i < problem.number_variables; i++) {
      if (is_finite_lower_bound(problem.variables_bounds[i].lb)) {
         error += std::abs(lower_bounds_multipliers[i] * (iterate.x[i] - problem.variables_bounds[i].lb));
      }
      if (is_finite_upper_bound(problem.variables_bounds[i].ub)) {
         error += std::abs(upper_bounds_multipliers[i] * (iterate.x[i] - problem.variables_bounds[i].ub));
      }
   }
   // constraints
   iterate.evaluate_constraints(problem, scaling);
   for (size_t j = 0; j < problem.number_constraints; j++) {
      const double multiplier_j = constraint_multipliers[j];
      const double scaled_lower_bound = scaling.get_constraint_scaling(j)*problem.constraint_bounds[j].lb;
      const double scaled_upper_bound = scaling.get_constraint_scaling(j)*problem.constraint_bounds[j].ub;
      if (iterate.constraints[j] < scaled_lower_bound) {
         // violated lower: the multiplier is 1 at optimum
         error += std::abs((1. - multiplier_j) * (scaled_lower_bound - iterate.constraints[j]));
      }
      else if (scaled_upper_bound < iterate.constraints[j]) {
         // violated upper: the multiplier is -1 at optimum
         error += std::abs((1. + multiplier_j) * (iterate.constraints[j] - scaled_upper_bound));
      }
      else if (is_finite_lower_bound(scaled_lower_bound) && 0. < multiplier_j) {
         error += std::abs(multiplier_j * (iterate.constraints[j] - scaled_lower_bound));
      }
      else if (is_finite_upper_bound(scaled_upper_bound) && multiplier_j < 0.) {
         error += std::abs(multiplier_j * (iterate.constraints[j] - scaled_upper_bound));
      }
   }
   return error;
}

void Subproblem::compute_optimality_conditions(const Problem& problem, const Scaling& scaling, Iterate& iterate, double objective_multiplier) const {
   iterate.evaluate_objective(problem, scaling);
   iterate.evaluate_constraints(problem, scaling);
   iterate.errors.constraints = problem.compute_constraint_violation(scaling, iterate.constraints, this->residual_norm);
   // compute the KKT error only if the objective multiplier is positive
   iterate.errors.KKT = this->compute_first_order_error(problem, scaling, iterate, 0. < objective_multiplier ? objective_multiplier : 1.);
   iterate.errors.FJ = this->compute_first_order_error(problem, scaling, iterate, 0.);
   iterate.errors.complementarity = Subproblem::compute_complementarity_error(problem, scaling, iterate, iterate.multipliers.constraints,
         iterate.multipliers.lower_bounds, iterate.multipliers.upper_bounds);
}

Direction Subproblem::compute_second_order_correction(const Problem& /*problem*/, Iterate& /*trial_iterate*/) {
   assert(false && "Subproblem::compute_second_order_correction");
}

void Subproblem::register_accepted_iterate(Iterate& /*iterate*/) {
   // by default, do nothing
}