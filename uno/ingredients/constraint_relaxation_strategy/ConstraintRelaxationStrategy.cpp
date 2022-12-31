// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"

ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Model& model, const Options& options):
      original_model(model),
      residual_norm(norm_from_string(options.get_string("residual_norm"))),
      small_step_threshold(options.get_double("small_step_threshold")) {
}

bool ConstraintRelaxationStrategy::is_small_step(const Direction& direction) const {
   // TODO for consistency with Uno.cpp, use tolerance/small_step_factor
   return (direction.norm <= this->small_step_threshold);
}

void ConstraintRelaxationStrategy::evaluate_lagrangian_gradient(Iterate& iterate, const std::vector<double>& constraint_multipliers,
      const std::vector<double>& lower_bound_multipliers, const std::vector<double>& upper_bound_multipliers) {
   initialize_vector(iterate.lagrangian_gradient, 0.);

   // objective gradient
   iterate.reformulation_evaluations.objective_gradient.for_each([&](size_t i, double derivative) {
         iterate.lagrangian_gradient[i] += derivative;
      });

   // constraints
   for (size_t j: Range(iterate.number_constraints)) {
      if (constraint_multipliers[j] != 0.) {
         iterate.reformulation_evaluations.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
            iterate.lagrangian_gradient[i] -= constraint_multipliers[j] * derivative;
         });
      }
   }

   // bound constraints
   for (size_t i: Range(iterate.number_variables)) {
      iterate.lagrangian_gradient[i] -= lower_bound_multipliers[i] + upper_bound_multipliers[i];
   }
}

void ConstraintRelaxationStrategy::evaluate_reformulation_functions(const NonlinearProblem& problem, Iterate& iterate) {
   // evaluate functions of the reformulated problem
   problem.evaluate_objective_gradient(iterate, iterate.reformulation_evaluations.objective_gradient);
   problem.evaluate_constraints(iterate, iterate.reformulation_evaluations.constraints);
   problem.evaluate_constraint_jacobian(iterate, iterate.reformulation_evaluations.constraint_jacobian);
}

double compute_stationarity_scaling(const NonlinearProblem& problem, const Iterate& iterate, double threshold) {
   const size_t total_size = problem.lower_bounded_variables.size() + problem.upper_bounded_variables.size() + problem.number_constraints;
   if (total_size == 0) {
      return 1.;
   }
   else {
      const double scaling_factor = threshold * static_cast<double>(total_size);
      return std::max(1., iterate.multipliers.compute_norm_1() / scaling_factor);
   }
}

double compute_complementarity_scaling(const NonlinearProblem& problem, const Iterate& iterate, double threshold) {
   const size_t total_size = problem.lower_bounded_variables.size() + problem.upper_bounded_variables.size();
   if (total_size == 0) {
      return 1.;
   }
   else {
      const double bound_multipliers_norm = norm_1(iterate.multipliers.lower_bounds) + norm_1(iterate.multipliers.upper_bounds);
      return std::max(1., bound_multipliers_norm / (threshold * static_cast<double>(total_size)));
   }
}

void ConstraintRelaxationStrategy::compute_primal_dual_errors(const NonlinearProblem& problem, Iterate& iterate) const {
   // stationarity error of the reformulated problem
   ConstraintRelaxationStrategy::evaluate_lagrangian_gradient(iterate, iterate.multipliers.constraints, iterate.multipliers.lower_bounds,
         iterate.multipliers.upper_bounds);
   iterate.residuals.stationarity = norm(iterate.lagrangian_gradient, this->residual_norm);
   // constraint violation of the original problem
   iterate.residuals.infeasibility = problem.model.compute_constraint_violation(iterate.model_evaluations.constraints, this->residual_norm);
   // complementarity error of the reformulated problem
   iterate.residuals.complementarity = problem.compute_complementarity_error(iterate.primals, iterate.reformulation_evaluations.constraints,
         iterate.multipliers.constraints, iterate.multipliers.lower_bounds, iterate.multipliers.upper_bounds);
   // scaling factors
   iterate.residuals.stationarity_scaling = compute_stationarity_scaling(problem, iterate, 100.);
   iterate.residuals.complementarity_scaling = compute_complementarity_scaling(problem, iterate, 100.);
   // TODO dual constraint violation of the reformulated problem
}

double ConstraintRelaxationStrategy::compute_linearized_constraint_violation(const Model& model, const Iterate& current_iterate,
      const Direction& direction, double step_length) {
   // determine the linearized constraint violation term: c(x_k) + alpha*\nabla c(x_k)^T d
   const auto jth_component = [&](size_t j) {
      const double component_j = current_iterate.model_evaluations.constraints[j] + step_length*dot(direction.primals,
            current_iterate.model_evaluations.constraint_jacobian[j]);
      return model.compute_constraint_violation(component_j, j);
   };

   return norm_1<double>(jth_component, Range(model.number_constraints));
}