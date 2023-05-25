// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"
#include "linear_algebra/view.hpp"

ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Model& model, const Options& options):
      original_model(model),
      progress_norm(norm_from_string(options.get_string("progress_norm"))),
      residual_norm(norm_from_string(options.get_string("residual_norm"))),
      residual_scaling_threshold(options.get_double("residual_scaling_threshold")) {
}

void ConstraintRelaxationStrategy::compute_primal_dual_residuals(const Model& model, const RelaxedProblem& feasibility_problem, Iterate& iterate) {
   iterate.evaluate_objective_gradient(model);
   iterate.evaluate_constraints(model);
   iterate.evaluate_constraint_jacobian(model);

   // stationarity error
   ConstraintRelaxationStrategy::evaluate_lagrangian_gradient(model.number_variables, iterate, iterate.multipliers, iterate.multipliers.objective);
   iterate.residuals.optimality_stationarity = this->compute_stationarity_error(iterate);
   iterate.residuals.feasibility_stationarity = feasibility_problem.compute_stationarity_error(iterate, this->residual_norm);

   // constraint violation of the original problem
   iterate.residuals.infeasibility = model.compute_constraint_violation(iterate.evaluations.constraints, this->residual_norm);

   // complementarity error
   iterate.residuals.optimality_complementarity = this->compute_complementarity_error(iterate.primals, iterate.evaluations.constraints,
         iterate.multipliers);
   iterate.residuals.feasibility_complementarity = feasibility_problem.compute_complementarity_error(iterate.primals, iterate.evaluations.constraints,
         iterate.multipliers, this->residual_norm);

   // scaling factors
   iterate.residuals.stationarity_scaling = this->compute_stationarity_scaling(model, iterate);
   iterate.residuals.complementarity_scaling = this->compute_complementarity_scaling(model, iterate);
}

// Lagrangian gradient split in two parts: objective contribution and constraints' contribution
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

double ConstraintRelaxationStrategy::compute_stationarity_error(const Iterate& iterate) const {
   // norm of the Lagrangian gradient
   return norm(this->residual_norm, iterate.lagrangian_gradient);
}

double ConstraintRelaxationStrategy::compute_stationarity_scaling(const Model& model, const Iterate& iterate) const {
   const size_t total_size = model.lower_bounded_variables.size() + model.upper_bounded_variables.size() + model.number_constraints;
   if (total_size == 0) {
      return 1.;
   }
   else {
      const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
      const double multiplier_norm = norm_1(
            view(iterate.multipliers.constraints, model.number_constraints),
            view(iterate.multipliers.lower_bounds, model.number_variables),
            view(iterate.multipliers.upper_bounds, model.number_variables)
      );
      return std::max(1., multiplier_norm / scaling_factor);
   }
}

double ConstraintRelaxationStrategy::compute_complementarity_scaling(const Model& model, const Iterate& iterate) const {
   const size_t total_size = model.lower_bounded_variables.size() + model.upper_bounded_variables.size();
   if (total_size == 0) {
      return 1.;
   }
   else {
      const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
      const double bound_multiplier_norm = norm_1(
            view(iterate.multipliers.lower_bounds, model.number_variables),
            view(iterate.multipliers.upper_bounds, model.number_variables)
      );
      return std::max(1., bound_multiplier_norm / scaling_factor);
   }
}