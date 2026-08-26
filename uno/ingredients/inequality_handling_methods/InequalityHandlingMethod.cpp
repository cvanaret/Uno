// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityHandlingMethod.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   InequalityHandlingMethod::InequalityHandlingMethod(const OptimizationProblem& problem, const Options& options):
         problem(problem), progress_norm(norm_from_string(options.get_string("progress_norm"))) {
   }

   // protected member functions

   void InequalityHandlingMethod::evaluate_progress_measures(const OptimizationProblem& problem, Iterate& iterate,
         Evaluations& evaluations) const {
      problem.set_infeasibility_measure(iterate, evaluations, this->progress_norm);
      problem.set_objective_measure(iterate, evaluations);
      problem.set_auxiliary_measure(iterate);
   }

   bool InequalityHandlingMethod::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Subproblem& subproblem, const SolverWorkspace& solver_workspace, const Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double step_length, Evaluations& current_evaluations, Evaluations& trial_evaluations) const {
      subproblem.problem.postprocess_iterate(trial_iterate);
      const double objective_multiplier = subproblem.problem.get_objective_multiplier();

      // evaluate progress measures
      evaluate_progress_measures(subproblem.problem, trial_iterate, trial_evaluations);
      trial_iterate.objective_multiplier = objective_multiplier;

      bool accept_iterate = false;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_evaluations.objective = this->problem.model.evaluate_objective(trial_iterate.primals);
         accept_iterate = true;
         statistics.set("Status", "0 primal step");
      }
      else {
         // determine acceptance wrt the globalization strategy
         const ProgressMeasures predicted_reductions = subproblem.compute_predicted_reductions(current_iterate, direction,
            step_length, this->progress_norm, current_evaluations, solver_workspace);
         accept_iterate = globalization_strategy.is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
            predicted_reductions, objective_multiplier);
         // check that the derivatives exist at the accepted trial iterate (an exception is thrown upon evaluation failure)
         if (accept_iterate) {
            trial_evaluations.evaluate_objective_gradient(this->problem.model, trial_iterate.primals);
            trial_evaluations.evaluate_jacobian(this->problem.model, trial_iterate.primals);
         }
      }
      return accept_iterate;
   }

   // stationarity errors:
   // - for KKT conditions: with standard multipliers and current objective multiplier
   // - for FJ conditions: with standard multipliers and 0 objective multiplier
   // - for feasibility problem: with feasibility multipliers and 0 objective multiplier
   void InequalityHandlingMethod::compute_residuals(const OptimizationProblem& problem, Iterate& iterate,
         Evaluations& evaluations) const {
      // stationarity error (norm of the Lagrangian gradient)
      problem.evaluate_lagrangian_gradient(iterate, evaluations, iterate.residuals.lagrangian_gradient);
      iterate.residuals.stationarity = norm(this->residual_norm, iterate.residuals.lagrangian_gradient);

      // primal feasibility/constraint violation of the model
      Vector<double> constraints(problem.number_constraints);
      problem.evaluate_constraints(iterate, constraints.data(), evaluations);
      iterate.primal_feasibility = problem.model.constraint_violation(evaluations.constraints, this->residual_norm);

      // complementarity error
      iterate.residuals.complementarity = problem.complementarity_error(iterate.primals, evaluations.constraints,
         iterate.multipliers, this->residual_norm);

      // scaling factors
      iterate.residuals.stationarity_scaling = this->compute_stationarity_scaling(problem.model, iterate.multipliers);
      iterate.residuals.complementarity_scaling = this->compute_complementarity_scaling(problem.model, iterate.multipliers);
   }

   double InequalityHandlingMethod::compute_stationarity_scaling(const Model& model, const Multipliers& multipliers) const {
      size_t number_lower_bounded_variables = 0;
      size_t number_upper_bounded_variables = 0;
      const auto& variables_lower_bounds = model.get_variables_lower_bounds();
      const auto& variables_upper_bounds = model.get_variables_upper_bounds();
      for (size_t variable_index: Range(model.number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index])) {
            ++number_lower_bounded_variables;
         }
         if (is_finite(variables_upper_bounds[variable_index])) {
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

   double InequalityHandlingMethod::compute_complementarity_scaling(const Model& model, const Multipliers& multipliers) const {
      size_t number_lower_bounded_variables = 0;
      size_t number_upper_bounded_variables = 0;
      const auto& variables_lower_bounds = model.get_variables_lower_bounds();
      const auto& variables_upper_bounds = model.get_variables_upper_bounds();
      for (size_t variable_index: Range(model.number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index])) {
            ++number_lower_bounded_variables;
         }
         if (is_finite(variables_upper_bounds[variable_index])) {
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