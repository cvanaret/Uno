// Copyright (c) 2018-2023 Charlie Vanaret
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

void ConstraintRelaxationStrategy::evaluate_lagrangian_gradient(size_t number_variables, Iterate& iterate, const Multipliers& multipliers,
      double objective_multiplier) {
   initialize_vector(iterate.lagrangian_gradient.objective_contribution, 0.);
   initialize_vector(iterate.lagrangian_gradient.constraints_contribution, 0.);

   // objective gradient
   iterate.evaluations.objective_gradient.for_each([&](size_t i, double derivative) {
         iterate.lagrangian_gradient.objective_contribution[i] += objective_multiplier * derivative;
      });

   // constraints
   for (size_t j: Range(iterate.number_constraints)) {
      if (multipliers.constraints[j] != 0.) {
         iterate.evaluations.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
            iterate.lagrangian_gradient.constraints_contribution[i] -= multipliers.constraints[j] * derivative;
         });
      }
   }

   // bound constraints
   for (size_t i: Range(number_variables)) {
      iterate.lagrangian_gradient.constraints_contribution[i] -= multipliers.lower_bounds[i] + multipliers.upper_bounds[i];
   }
}

double compute_stationarity_scaling(const NonlinearProblem& problem, const Iterate& iterate, double threshold) {
   const size_t total_size = problem.lower_bounded_variables.size() + problem.upper_bounded_variables.size() + problem.number_constraints;
   if (total_size == 0) {
      return 1.;
   }
   else {
      const double scaling_factor = threshold * static_cast<double>(total_size);
      return std::max(1., iterate.multipliers.norm_1() / scaling_factor);
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

void ConstraintRelaxationStrategy::compute_primal_dual_residuals(const NonlinearProblem& problem, Iterate& iterate) const {
   iterate.evaluate_objective_gradient(problem.model);
   iterate.evaluate_constraints(problem.model);
   iterate.evaluate_constraint_jacobian(problem.model);

   // stationarity error
   ConstraintRelaxationStrategy::evaluate_lagrangian_gradient(problem.model.number_variables, iterate, iterate.multipliers,
         problem.get_objective_multiplier());
   iterate.residuals.optimality_stationarity = NonlinearProblem::compute_optimality_stationarity_error(iterate, this->residual_norm);
   iterate.residuals.feasibility_stationarity = NonlinearProblem::compute_feasibility_stationarity_error(iterate, this->residual_norm);

   // constraint violation of the original problem
   iterate.residuals.infeasibility = problem.model.compute_constraint_violation(iterate.evaluations.constraints, this->residual_norm);

   // complementarity error
   iterate.residuals.optimality_complementarity = problem.compute_complementarity_error(problem.model.number_variables, iterate.primals,
         iterate.evaluations.constraints, iterate.multipliers);
   iterate.residuals.feasibility_complementarity = problem.compute_feasibility_complementarity_error(problem.model.number_variables, iterate.primals,
         iterate.evaluations.constraints, iterate.multipliers);

   // scaling factors
   // TODO put these coefficients in the option file
   iterate.residuals.stationarity_scaling = compute_stationarity_scaling(problem, iterate, 100.);
   iterate.residuals.complementarity_scaling = compute_complementarity_scaling(problem, iterate, 100.);
   // TODO dual constraint violation of the reformulated problem
}

double ConstraintRelaxationStrategy::compute_linearized_constraint_violation(const Model& model, const Iterate& current_iterate,
      const Direction& direction, double step_length) {
   // determine the linearized constraint violation term: c(x_k) + alpha*\nabla c(x_k)^T d
   const auto jth_component = [&](size_t j) {
      const double component_j = current_iterate.evaluations.constraints[j] + step_length * dot(direction.primals,
            current_iterate.evaluations.constraint_jacobian[j]);
      return model.compute_constraint_violation(component_j, j);
   };

   return norm_1<double>(jth_component, Range(model.number_constraints));
}