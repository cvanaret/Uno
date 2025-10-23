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

namespace uno {
   ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Options& options):
         residual_norm(norm_from_string(options.get_string("residual_norm"))),
         residual_scaling_threshold(options.get_double("residual_scaling_threshold")),
         primal_tolerance(options.get_double("primal_tolerance")),
         dual_tolerance(options.get_double("dual_tolerance")),
         loose_primal_tolerance(options.get_double("loose_primal_tolerance")),
         loose_dual_tolerance(options.get_double("loose_dual_tolerance")),
         loose_tolerance_consecutive_iteration_threshold(options.get_unsigned_int("loose_tolerance_consecutive_iteration_threshold")),
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