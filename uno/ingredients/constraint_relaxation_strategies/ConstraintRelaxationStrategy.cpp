// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/VectorView.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Multipliers.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/Options.hpp"
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
         loose_primal_tolerance(options.get_double("loose_primal_tolerance")),
         loose_dual_tolerance(options.get_double("loose_dual_tolerance")),
         loose_tolerance_iteration_threshold(options.get_unsigned_int("loose_tolerance_iteration_threshold")),
         unbounded_objective_threshold(options.get_double("unbounded_objective_threshold")) {
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

   void ConstraintRelaxationStrategy::evaluate_progress_measures(const OptimizationProblem& problem, Iterate& iterate,
         Evaluations& evaluations) const {
      problem.set_infeasibility_measure(iterate, evaluations, this->progress_norm);
      problem.set_objective_measure(iterate, evaluations);
      problem.set_auxiliary_measure(iterate);
   }

   bool ConstraintRelaxationStrategy::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Subproblem& subproblem, const SolverWorkspace& solver_workspace, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double step_length, EvaluationCache& evaluation_cache, UserCallbacks& user_callbacks) const {
      subproblem.problem.postprocess_iterate(trial_iterate);
      const double objective_multiplier = subproblem.problem.get_objective_multiplier();

      // evaluate progress measures
      trial_iterate.objective_multiplier = objective_multiplier;
      //if (this->subproblem_definition_changed) {
         //DEBUG << "The subproblem definition changed, the globalization strategy is reset and the auxiliary measure is recomputed\n";
         //globalization_strategy.reset();
         subproblem.problem.set_auxiliary_measure(current_iterate);
         //this->subproblem_definition_changed = false;
      //}
      this->evaluate_progress_measures(subproblem.problem, trial_iterate, evaluation_cache.trial_evaluations);

      bool accept_iterate = false;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         evaluation_cache.trial_evaluations.objective = subproblem.problem.model.evaluate_objective(trial_iterate.primals);
         accept_iterate = true;
         statistics.set("Status", "0 primal step");
      }
      else {
         // determine acceptance wrt the globalization strategy
         const ProgressMeasures predicted_reductions = subproblem.compute_predicted_reductions(direction, step_length,
            this->progress_norm, evaluation_cache.current_evaluations, solver_workspace);
         accept_iterate = globalization_strategy.is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
            predicted_reductions, objective_multiplier);
         // check that the derivatives exist at the accepted trial iterate (an exception is thrown upon evaluation failure)
         if (accept_iterate) {
            evaluation_cache.trial_evaluations.evaluate_objective_gradient(subproblem.problem.model, trial_iterate.primals);
            evaluation_cache.trial_evaluations.evaluate_jacobian(subproblem.problem.model, trial_iterate.primals);
         }
      }
      if (accept_iterate) {
         user_callbacks.notify_acceptable_iterate(trial_iterate.primals, trial_iterate.multipliers, objective_multiplier,
            trial_iterate.progress.infeasibility, trial_iterate.residuals.stationarity, trial_iterate.residuals.complementarity);
      }
      return accept_iterate;
   }

   // stationarity errors:
   // - for KKT conditions: with standard multipliers and current objective multiplier
   // - for FJ conditions: with standard multipliers and 0 objective multiplier
   // - for feasibility problem: with feasibility multipliers and 0 objective multiplier
   void ConstraintRelaxationStrategy::compute_residuals(const OptimizationProblem& problem, Iterate& iterate,
         Evaluations& evaluations) const {
      // stationarity error (norm of the Lagrangian gradient)
      problem.evaluate_lagrangian_gradient(iterate, evaluations, iterate.residuals.lagrangian_gradient);
      iterate.residuals.stationarity = norm(this->residual_norm, iterate.residuals.lagrangian_gradient);

      // primal feasibility/constraint violation of the model
      evaluations.evaluate_constraints(problem.model, iterate.primals);
      iterate.primal_feasibility = problem.model.constraint_violation(evaluations.constraints, this->residual_norm);

      // complementarity error
      constexpr double shift_value = 0.;
      // TODO preallocate constraints
      Vector<double> constraints(problem.number_constraints);
      problem.evaluate_constraints(iterate, constraints.data(), evaluations);
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