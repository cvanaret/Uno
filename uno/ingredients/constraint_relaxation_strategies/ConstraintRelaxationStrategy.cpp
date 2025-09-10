// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Multipliers.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"
#include "symbolic/Expression.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Options& options):
         progress_norm(norm_from_string(options.get_string("progress_norm"))),
         residual_norm(norm_from_string(options.get_string("residual_norm"))),
         residual_scaling_threshold(options.get_double("residual_scaling_threshold")),
         primal_tolerance(options.get_double("primal_tolerance")),
         dual_tolerance(options.get_double("dual_tolerance")),
         loose_dual_tolerance(options.get_double("loose_dual_tolerance")),
         loose_tolerance_consecutive_iteration_threshold(options.get_unsigned_int("loose_tolerance_consecutive_iteration_threshold")),
         unbounded_objective_threshold(options.get_double("unbounded_objective_threshold")),
         first_order_predicted_reduction(options.get_string("globalization_mechanism") == "LS") {
   }

   ConstraintRelaxationStrategy::~ConstraintRelaxationStrategy() { }

   // with initial point
   /*
   void ConstraintRelaxationStrategy::compute_feasible_direction(Statistics& statistics, InequalityHandlingMethod& inequality_handling_method,
         GlobalizationStrategy& globalization_strategy, const Model& model, Iterate& current_iterate, Direction& direction,
         const Vector<double>& initial_point, double trust_region_radius, WarmstartInformation& warmstart_information) {
      inequality_handling_method.set_initial_point(initial_point);
      this->compute_feasible_direction(statistics, globalization_strategy, model, current_iterate,
         direction, trust_region_radius, warmstart_information);
   }
   */

   // infeasibility measure: constraint violation
   void ConstraintRelaxationStrategy::set_infeasibility_measure(const Model& model, Iterate& iterate) const {
      iterate.evaluate_constraints(model);
      iterate.progress.infeasibility = model.constraint_violation(iterate.evaluations.constraints, this->progress_norm);
   }

   // objective measure: scaled objective
   void ConstraintRelaxationStrategy::set_objective_measure(const Model& model, Iterate& iterate) const {
      iterate.evaluate_objective(model);
      const double objective = iterate.evaluations.objective;
      iterate.progress.objective = [=](double objective_multiplier) {
         return objective_multiplier * objective;
      };
   }

   double ConstraintRelaxationStrategy::compute_predicted_infeasibility_reduction(InequalityHandlingMethod& inequality_handling_method,
         const Model& model, const Iterate& current_iterate, const Vector<double>& primal_direction, double step_length) const {
      // predicted infeasibility reduction: "‖c(x)‖ - ‖c(x) + ∇c(x)^T (αd)‖"
      const double current_constraint_violation = model.constraint_violation(current_iterate.evaluations.constraints, this->progress_norm);
      Vector<double> result(model.number_constraints);
      inequality_handling_method.compute_constraint_jacobian_vector_product(primal_direction, result);
      const double trial_linearized_constraint_violation = model.constraint_violation(current_iterate.evaluations.constraints +
         step_length * result, this->progress_norm);
      return current_constraint_violation - trial_linearized_constraint_violation;
   }

   std::function<double(double)> ConstraintRelaxationStrategy::compute_predicted_objective_reduction(InequalityHandlingMethod& inequality_handling_method,
         const Iterate& current_iterate, const Vector<double>& primal_direction, double step_length) const {
      // predicted objective reduction: "-∇f(x)^T (αd) - α^2/2 d^T H d"
      const double directional_derivative = dot(primal_direction, current_iterate.evaluations.objective_gradient);
      const double quadratic_term = this->first_order_predicted_reduction ? 0. :
         inequality_handling_method.compute_hessian_quadratic_product(primal_direction);
      return [=](double objective_multiplier) {
         return step_length * (-objective_multiplier*directional_derivative) - step_length*step_length/2. * quadratic_term;
      };
   }

   void ConstraintRelaxationStrategy::compute_progress_measures(InequalityHandlingMethod& inequality_handling_method,
         const OptimizationProblem& problem, GlobalizationStrategy& globalization_strategy, Iterate& current_iterate,
         Iterate& trial_iterate) const {
      if (inequality_handling_method.subproblem_definition_changed) {
         DEBUG << "The subproblem definition changed, the globalization strategy is reset and the auxiliary measure is recomputed\n";
         globalization_strategy.reset();
         inequality_handling_method.set_auxiliary_measure(problem, current_iterate);
         inequality_handling_method.subproblem_definition_changed = false;
      }
      this->evaluate_progress_measures(inequality_handling_method, problem, trial_iterate);
   }

   ProgressMeasures ConstraintRelaxationStrategy::compute_predicted_reductions(InequalityHandlingMethod& inequality_handling_method,
         const OptimizationProblem& problem, const Iterate& current_iterate, const Direction& direction, double step_length) const {
      return {
         this->compute_predicted_infeasibility_reduction(inequality_handling_method, problem.model, current_iterate, direction.primals, step_length),
         this->compute_predicted_objective_reduction(inequality_handling_method, current_iterate, direction.primals, step_length),
         inequality_handling_method.compute_predicted_auxiliary_reduction_model(problem, current_iterate, direction.primals, step_length)
      };
   }

   bool ConstraintRelaxationStrategy::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const OptimizationProblem& problem, InequalityHandlingMethod& inequality_handling_method, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, UserCallbacks& user_callbacks) const {
      inequality_handling_method.postprocess_iterate(problem, trial_iterate);
      const double objective_multiplier = problem.get_objective_multiplier();
      trial_iterate.objective_multiplier = objective_multiplier;
      this->compute_progress_measures(inequality_handling_method, problem, globalization_strategy, current_iterate, trial_iterate);

      bool accept_iterate = false;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_iterate.evaluate_objective(problem.model);
         accept_iterate = true;
         statistics.set("status", "0 primal step");
      }
      else {
         const ProgressMeasures predicted_reductions = ConstraintRelaxationStrategy::compute_predicted_reductions(inequality_handling_method,
            problem, current_iterate, direction, step_length);
         accept_iterate = globalization_strategy.is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
            predicted_reductions, objective_multiplier);
      }
      if (accept_iterate) {
         user_callbacks.notify_acceptable_iterate(trial_iterate.primals, trial_iterate.multipliers, objective_multiplier);
      }
      return accept_iterate;
   }

   // stationarity errors:
   // - for KKT conditions: with standard multipliers and current objective multiplier
   // - for FJ conditions: with standard multipliers and 0 objective multiplier
   // - for feasibility problem: with feasibility multipliers and 0 objective multiplier

   void ConstraintRelaxationStrategy::compute_primal_dual_residuals(const OptimizationProblem& problem, Iterate& iterate) const {
      iterate.residuals.stationarity = OptimizationProblem::stationarity_error(iterate.residuals.lagrangian_gradient,
         problem.get_objective_multiplier(), this->residual_norm);

      // constraint violation of the original problem
      iterate.primal_feasibility = problem.model.constraint_violation(iterate.evaluations.constraints, this->residual_norm);

      // complementarity error
      constexpr double shift_value = 0.;
      // TODO preallocate constraints
      Vector<double> constraints(problem.number_constraints);
      problem.evaluate_constraints(iterate, constraints);
      iterate.residuals.complementarity = problem.complementarity_error(iterate.primals, constraints,
         iterate.multipliers, shift_value, this->residual_norm);

      // scaling factors
      iterate.residuals.stationarity_scaling = this->compute_stationarity_scaling(problem.model, iterate.multipliers);
      iterate.residuals.complementarity_scaling = this->compute_complementarity_scaling(problem.model, iterate.multipliers);
   }

   double ConstraintRelaxationStrategy::compute_stationarity_scaling(const Model& model, const Multipliers& multipliers) const {
      size_t number_lower_bounded_variables = 0;
      size_t number_upper_bounded_variables = 0;
      for (size_t variable_index: Range(model.number_variables)) {
         if (is_finite(model.variable_lower_bound(variable_index))) {
            ++number_lower_bounded_variables;
         }
         if (is_finite(model.variable_upper_bound(variable_index))) {
            ++number_upper_bounded_variables;
         }
      }
      const size_t total_size = number_lower_bounded_variables + number_upper_bounded_variables + model.number_constraints;
      if (total_size == 0) {
         return 1.;
      }
      else {
         const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
         const double multiplier_norm = norm_1(
               view(multipliers.constraints, 0, model.number_constraints),
               view(multipliers.lower_bounds, 0, model.number_variables),
               view(multipliers.upper_bounds, 0, model.number_variables)
         );
         return std::max(1., multiplier_norm / scaling_factor);
      }
   }

   double ConstraintRelaxationStrategy::compute_complementarity_scaling(const Model& model, const Multipliers& multipliers) const {
      size_t number_lower_bounded_variables = 0;
      size_t number_upper_bounded_variables = 0;
      for (size_t variable_index: Range(model.number_variables)) {
         if (is_finite(model.variable_lower_bound(variable_index))) {
            ++number_lower_bounded_variables;
         }
         if (is_finite(model.variable_upper_bound(variable_index))) {
            ++number_upper_bounded_variables;
         }
      }
      const size_t total_size = number_lower_bounded_variables + number_upper_bounded_variables;
      if (total_size == 0) {
         return 1.;
      }
      else {
         const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
         const double bound_multiplier_norm = norm_1(
               view(multipliers.lower_bounds, 0, model.number_variables),
               view(multipliers.upper_bounds, 0, model.number_variables)
         );
         return std::max(1., bound_multiplier_norm / scaling_factor);
      }
   }
} // namespace