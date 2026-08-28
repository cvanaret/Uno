// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/View.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Multipliers.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Options& options):
         residual_norm(norm_from_string(options.get_string("residual_norm"))),
         residual_scaling_threshold(options.get_double("residual_scaling_threshold")),
         primal_tolerance(options.get_double("primal_tolerance")),
         dual_tolerance(options.get_double("dual_tolerance")),
         loose_primal_tolerance(options.get_double("loose_primal_tolerance")),
         loose_dual_tolerance(options.get_double("loose_dual_tolerance")),
         loose_tolerance_iteration_threshold(options.get_unsigned_int("loose_tolerance_iteration_threshold")),
         diverging_iterate_threshold(options.get_double("diverging_iterate_threshold")),
         unbounded_objective_threshold(options.get_double("unbounded_objective_threshold")) {
   }

   ConstraintRelaxationStrategy::~ConstraintRelaxationStrategy() = default;

   size_t ConstraintRelaxationStrategy::get_number_subproblems_solved() const {
      return this->number_subproblems_solved;
   }

   // protected member functions

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
      problem.evaluate_constraints(iterate, constraints.view(), evaluations);
      iterate.residuals.complementarity = problem.complementarity_error(iterate.primals, constraints,
         iterate.multipliers, shift_value, this->residual_norm);

      // scaling factors
      iterate.residuals.stationarity_scaling = this->compute_stationarity_scaling(problem.model, iterate.multipliers);
      iterate.residuals.complementarity_scaling = this->compute_complementarity_scaling(problem.model, iterate.multipliers);
   }

   double ConstraintRelaxationStrategy::compute_stationarity_scaling(const Model& model, const Multipliers& multipliers) const {
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

   double ConstraintRelaxationStrategy::compute_complementarity_scaling(const Model& model, const Multipliers& multipliers) const {
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