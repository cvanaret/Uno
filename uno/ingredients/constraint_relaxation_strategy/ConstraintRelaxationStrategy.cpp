// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/LagrangianGradient.hpp"
#include "optimization/Multipliers.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "symbolic/VectorView.hpp"
#include "symbolic/Expression.hpp"
#include "tools/Options.hpp"
#include "tools/Statistics.hpp"

ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Model& model, size_t number_variables, size_t number_constraints,
      size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros, const Options& options):
      model(model),
      globalization_strategy(GlobalizationStrategyFactory::create(options.get_string("globalization_strategy"), options)),
      subproblem(SubproblemFactory::create(number_variables, number_constraints, number_objective_gradient_nonzeros, number_jacobian_nonzeros,
            number_hessian_nonzeros, options)),
      progress_norm(norm_from_string(options.get_string("progress_norm"))),
      residual_norm(norm_from_string(options.get_string("residual_norm"))),
      residual_scaling_threshold(options.get_double("residual_scaling_threshold")),
      tight_tolerance(options.get_double("tolerance")),
      loose_tolerance(options.get_double("loose_tolerance")),
      loose_tolerance_consecutive_iteration_threshold(options.get_unsigned_int("loose_tolerance_consecutive_iteration_threshold")),
      unbounded_objective_threshold(options.get_double("unbounded_objective_threshold")) {
}

void ConstraintRelaxationStrategy::set_trust_region_radius(double trust_region_radius) {
   this->subproblem->set_trust_region_radius(trust_region_radius);
}

// with initial point
void ConstraintRelaxationStrategy::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
      const Vector<double>& initial_point, WarmstartInformation& warmstart_information) {
   this->subproblem->set_initial_point(initial_point);
   this->compute_feasible_direction(statistics, current_iterate, direction, warmstart_information);
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

double ConstraintRelaxationStrategy::compute_predicted_infeasibility_reduction_model(const Iterate& current_iterate,
      const Vector<double>& primal_direction, double step_length) const {
   // predicted infeasibility reduction: "‖c(x)‖ - ‖c(x) + ∇c(x)^T (αd)‖"
   const double current_constraint_violation = this->model.constraint_violation(current_iterate.evaluations.constraints, this->progress_norm);
   const double trial_linearized_constraint_violation = this->model.constraint_violation(current_iterate.evaluations.constraints + step_length *
         (current_iterate.evaluations.constraint_jacobian * primal_direction), this->progress_norm);
   return current_constraint_violation - trial_linearized_constraint_violation;
}

std::function<double(double)> ConstraintRelaxationStrategy::compute_predicted_objective_reduction_model(const Iterate& current_iterate,
      const Vector<double>& primal_direction, double step_length, const SymmetricMatrix<double>& hessian) const {
   // predicted objective reduction: "-∇f(x)^T (αd) - α^2/2 d^T H d"
   const double directional_derivative = dot(primal_direction, current_iterate.evaluations.objective_gradient);
   const double quadratic_term = hessian.quadratic_product(primal_direction, primal_direction);
   return [=](double objective_multiplier) {
      return step_length * (-objective_multiplier*directional_derivative) - step_length*step_length/2. * quadratic_term;
   };
}

void ConstraintRelaxationStrategy::compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate) {
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed, the globalization strategy is reset and the auxiliary measure is recomputed\n";
      this->globalization_strategy->reset();
      this->subproblem->set_auxiliary_measure(this->model, current_iterate);
      this->subproblem->subproblem_definition_changed = false;
   }
   this->evaluate_progress_measures(trial_iterate);
}

void ConstraintRelaxationStrategy::compute_primal_dual_residuals(const OptimizationProblem& problem, Iterate& iterate, const Multipliers& multipliers) {
   iterate.evaluate_objective_gradient(this->model);
   iterate.evaluate_constraints(this->model);
   iterate.evaluate_constraint_jacobian(this->model);

   // stationarity errors:
   // - for KKT conditions: with standard multipliers and current objective multiplier
   // - for FJ conditions: with standard multipliers and 0 objective multiplier
   // - for feasibility problem: with feasibility multipliers and 0 objective multiplier
   problem.evaluate_lagrangian_gradient(iterate, multipliers);
   iterate.residuals.stationarity = ConstraintRelaxationStrategy::stationarity_error(iterate.lagrangian_gradient, problem.get_objective_multiplier(),
         this->residual_norm);

   // constraint violation of the original problem
   iterate.residuals.primal_feasibility = this->model.constraint_violation(iterate.evaluations.constraints, this->residual_norm);

   iterate.residuals.dual_feasibility = problem.dual_feasibility_error(multipliers, this->residual_norm);

   // complementarity error
   const double shift_value = 0.;
   iterate.residuals.complementarity = problem.complementarity_error(iterate.primals, iterate.evaluations.constraints, multipliers, shift_value,
         this->residual_norm);

   // scaling factors
   iterate.residuals.stationarity_scaling = this->compute_stationarity_scaling(multipliers);
   iterate.residuals.complementarity_scaling = this->compute_complementarity_scaling(multipliers);
}

double ConstraintRelaxationStrategy::stationarity_error(const LagrangianGradient<double>& lagrangian_gradient, double objective_multiplier,
      Norm residual_norm) {
   // norm of the scaled Lagrangian gradient
   const auto scaled_lagrangian = objective_multiplier * lagrangian_gradient.objective_contribution + lagrangian_gradient.constraints_contribution;
   return norm(residual_norm, scaled_lagrangian);
}

double ConstraintRelaxationStrategy::compute_stationarity_scaling(const Multipliers& multipliers) const {
   const size_t total_size = this->model.get_lower_bounded_variables().size() + this->model.get_upper_bounded_variables().size() +
         this->model.number_constraints;
   if (total_size == 0) {
      return 1.;
   }
   else {
      const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
      const double multiplier_norm = norm_1(
            view(multipliers.constraints, 0, this->model.number_constraints),
            view(multipliers.lower_bounds, 0, this->model.number_variables),
            view(multipliers.upper_bounds, 0, this->model.number_variables)
      );
      return std::max(1., multiplier_norm / scaling_factor);
   }
}

double ConstraintRelaxationStrategy::compute_complementarity_scaling(const Multipliers& multipliers) const {
   const size_t total_size = this->model.get_lower_bounded_variables().size() + this->model.get_upper_bounded_variables().size();
   if (total_size == 0) {
      return 1.;
   }
   else {
      const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
      const double bound_multiplier_norm = norm_1(
            view(multipliers.lower_bounds, 0, this->model.number_variables),
            view(multipliers.upper_bounds, 0, this->model.number_variables)
      );
      return std::max(1., bound_multiplier_norm / scaling_factor);
   }
}

TerminationStatus ConstraintRelaxationStrategy::check_termination(const OptimizationProblem& problem, Iterate& current_iterate) {
   // test convergence wrt the tight tolerance
   const TerminationStatus status_tight_tolerance = problem.check_convergence_with_given_tolerance(current_iterate, this->tight_tolerance);
   if (status_tight_tolerance != TerminationStatus::NOT_OPTIMAL || this->loose_tolerance <= this->tight_tolerance) {
      return status_tight_tolerance;
   }

   // if not converged, check convergence wrt loose tolerance (provided it is strictly looser than the tight tolerance)
   const TerminationStatus status_loose_tolerance = problem.check_convergence_with_given_tolerance(current_iterate, this->loose_tolerance);
   // if converged, keep track of the number of consecutive iterations
   if (status_loose_tolerance != TerminationStatus::NOT_OPTIMAL) {
      this->loose_tolerance_consecutive_iterations++;
   }
   else {
      this->loose_tolerance_consecutive_iterations = 0;
      return TerminationStatus::NOT_OPTIMAL;
   }
   // check if loose tolerance achieved for enough consecutive iterations
   if (this->loose_tolerance_consecutive_iteration_threshold <= this->loose_tolerance_consecutive_iterations) {
      return status_loose_tolerance;
   }
   else {
      return TerminationStatus::NOT_OPTIMAL;
   }
}

void ConstraintRelaxationStrategy::set_statistics(Statistics& statistics, const Iterate& iterate) const {
   this->set_progress_statistics(statistics, iterate);
   this->set_dual_residuals_statistics(statistics, iterate);
}

void ConstraintRelaxationStrategy::set_progress_statistics(Statistics& statistics, const Iterate& iterate) const {
   statistics.set("objective", iterate.evaluations.objective);
   if (this->model.is_constrained()) {
      statistics.set("primal feas.", iterate.progress.infeasibility);
   }
}

size_t ConstraintRelaxationStrategy::get_hessian_evaluation_count() const {
   return this->subproblem->get_hessian_evaluation_count();
}

size_t ConstraintRelaxationStrategy::get_number_subproblems_solved() const {
   return this->subproblem->number_subproblems_solved;
}