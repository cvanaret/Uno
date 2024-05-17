// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Multipliers.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "symbolic/VectorView.hpp"
#include "symbolic/Expression.hpp"
#include "tools/Options.hpp"
#include "tools/Statistics.hpp"

ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Model& model, const Options& options):
      model(model),
      progress_norm(norm_from_string(options.get_string("progress_norm"))),
      residual_norm(norm_from_string(options.get_string("residual_norm"))),
      residual_scaling_threshold(options.get_double("residual_scaling_threshold")) {
}

// infeasibility measure: constraint violation
void ConstraintRelaxationStrategy::set_infeasibility_measure(Iterate& iterate) const {
   iterate.evaluate_constraints(this->model);
   iterate.progress.infeasibility = this->model.constraint_violation(iterate.evaluations.constraints, this->progress_norm);
}

// objective measure: scaled objective
void ConstraintRelaxationStrategy::set_objective_measure(Iterate& iterate) const {
   iterate.evaluate_objective(this->model);
   const double objective = iterate.evaluations.objective;
   iterate.progress.objective = [=](double objective_multiplier) {
      return objective_multiplier * objective;
   };
}

double ConstraintRelaxationStrategy::compute_predicted_infeasibility_reduction_model(const Iterate& current_iterate, const Direction& direction,
      double step_length) const {
   // predicted infeasibility reduction: "‖c(x)‖ - ‖c(x) + ∇c(x)^T (αd)‖"
   const double current_constraint_violation = this->model.constraint_violation(current_iterate.evaluations.constraints, this->progress_norm);
   const double trial_linearized_constraint_violation = this->model.constraint_violation(current_iterate.evaluations.constraints + step_length *
         (current_iterate.evaluations.constraint_jacobian * direction.primals), this->progress_norm);
   return current_constraint_violation - trial_linearized_constraint_violation;
}

std::function<double(double)> ConstraintRelaxationStrategy::compute_predicted_objective_reduction_model(const Iterate& current_iterate,
      const Direction& direction, double step_length, const SymmetricMatrix<double>& hessian) const {
   // predicted objective reduction: "-∇f(x)^T (αd) - α^2/2 d^T H d"
   const double directional_derivative = dot(direction.primals, current_iterate.evaluations.objective_gradient);
   const double quadratic_term = hessian.quadratic_product(direction.primals, direction.primals);
   return [=](double objective_multiplier) {
      return step_length * (-objective_multiplier*directional_derivative) - step_length*step_length/2. * quadratic_term;
   };
}

void ConstraintRelaxationStrategy::compute_primal_dual_residuals(const OptimizationProblem& optimality_problem, const OptimizationProblem& feasibility_problem,
      Iterate& iterate) {
   iterate.evaluate_objective_gradient(this->model);
   iterate.evaluate_constraints(this->model);
   iterate.evaluate_constraint_jacobian(this->model);

   // stationarity error
   this->evaluate_lagrangian_gradient(iterate, iterate.multipliers);
   iterate.residuals.optimality_stationarity = optimality_problem.stationarity_error(iterate.lagrangian_gradient, iterate.objective_multiplier,
         this->residual_norm);
   iterate.residuals.feasibility_stationarity = feasibility_problem.stationarity_error(iterate.lagrangian_gradient, 0., this->residual_norm);

   // constraint violation of the original problem
   iterate.residuals.infeasibility = this->model.constraint_violation(iterate.evaluations.constraints, this->residual_norm);

   // complementarity error
   iterate.residuals.optimality_complementarity = optimality_problem.complementarity_error(iterate.primals, iterate.evaluations.constraints,
         iterate.multipliers, this->residual_norm);
   iterate.residuals.feasibility_complementarity = feasibility_problem.complementarity_error(iterate.primals, iterate.evaluations.constraints,
         iterate.multipliers, this->residual_norm);

   // scaling factors
   iterate.residuals.stationarity_scaling = this->compute_stationarity_scaling(iterate);
   iterate.residuals.complementarity_scaling = this->compute_complementarity_scaling(iterate);
}

// Lagrangian gradient split in two parts: objective contribution and constraints' contribution
void ConstraintRelaxationStrategy::evaluate_lagrangian_gradient(Iterate& iterate, const Multipliers& multipliers) const {
   initialize_vector(iterate.lagrangian_gradient.objective_contribution, 0.);
   initialize_vector(iterate.lagrangian_gradient.constraints_contribution, 0.);

   // objective gradient
   for (auto [variable_index, derivative]: iterate.evaluations.objective_gradient) {
      iterate.lagrangian_gradient.objective_contribution[variable_index] += derivative;
   }

   // constraints
   for (size_t constraint_index: Range(iterate.number_constraints)) {
      if (multipliers.constraints[constraint_index] != 0.) {
         for (auto [variable_index, derivative]: iterate.evaluations.constraint_jacobian[constraint_index]) {
            iterate.lagrangian_gradient.constraints_contribution[variable_index] -= multipliers.constraints[constraint_index] * derivative;
         }
      }
   }

   // bound constraints
   for (size_t variable_index: Range(this->model.number_variables)) {
      iterate.lagrangian_gradient.constraints_contribution[variable_index] -= multipliers.lower_bounds[variable_index] + multipliers.upper_bounds[variable_index];
   }
}

double ConstraintRelaxationStrategy::compute_stationarity_scaling(const Iterate& iterate) const {
   const size_t total_size = this->model.get_lower_bounded_variables().size() + this->model.get_upper_bounded_variables().size() + this->model.number_constraints;
   if (total_size == 0) {
      return 1.;
   }
   else {
      const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
      const double multiplier_norm = norm_1(
            view(iterate.multipliers.constraints, this->model.number_constraints),
            view(iterate.multipliers.lower_bounds, this->model.number_variables),
            view(iterate.multipliers.upper_bounds, this->model.number_variables)
      );
      return std::max(1., multiplier_norm / scaling_factor);
   }
}

double ConstraintRelaxationStrategy::compute_complementarity_scaling(const Iterate& iterate) const {
   const size_t total_size = this->model.get_lower_bounded_variables().size() + this->model.get_upper_bounded_variables().size();
   if (total_size == 0) {
      return 1.;
   }
   else {
      const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
      const double bound_multiplier_norm = norm_1(
            view(iterate.multipliers.lower_bounds, this->model.number_variables),
            view(iterate.multipliers.upper_bounds, this->model.number_variables)
      );
      return std::max(1., bound_multiplier_norm / scaling_factor);
   }
}

void ConstraintRelaxationStrategy::set_statistics(Statistics& statistics, const Iterate& iterate) const {
   this->set_progress_statistics(statistics, iterate);
   this->set_dual_residuals_statistics(statistics, iterate);
}

void ConstraintRelaxationStrategy::set_progress_statistics(Statistics& statistics, const Iterate& iterate) const {
   statistics.set("objective", iterate.evaluations.objective);
   if (this->model.is_constrained()) {
      statistics.set("infeasibility", iterate.progress.infeasibility);
   }
}