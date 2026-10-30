// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cmath>
#include <exception>
#include <random>
#include <utility>
#include "ModelReduction.hpp"
#include "Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "options/Options.hpp"
#include "symbolic/Collection.hpp"
#include "symbolic/Range.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   namespace {
      // term J_e (x_i - x̄_i) of a linear constraint c_j(x) = c_j(x̄) + sum_e J_e (x_i - x̄_i). Duplicates are allowed
      struct LinearTerm {
         size_t variable_index;
         double coefficient;
      };
      using LinearRows = std::vector<std::vector<LinearTerm>>;

      bool are_close(double value1, double value2) {
         return std::abs(value1 - value2) <= 1e-12 * std::max(1., std::abs(value1));
      }

      // nonlinear constraints whose Jacobian at a couple of random points around the base point (within the bounds) is
      // that at the base point. The fixed variables stay at their values: the constraints are linear in the reduced space
      std::vector<size_t> detect_linear_constraints(const Model& model, const Vector<double>& base_point,
            const Vector<double>& base_jacobian, const std::vector<bool>& is_fixed) {
         constexpr size_t number_random_points = 2;
         std::vector<bool> is_constant(model.number_constraints, false);
         for (size_t constraint_index: model.get_nonlinear_constraints()) {
            is_constant[constraint_index] = true;
         }
         std::mt19937 generator(20261030); // fixed seed: reproducible
         std::uniform_real_distribution<double> distribution(-1., 1.);
         const auto& lower_bounds = model.get_variables_lower_bounds();
         const auto& upper_bounds = model.get_variables_upper_bounds();
         const auto row_indices = model.get_jacobian_row_indices();
         const auto column_indices = model.get_jacobian_column_indices();
         Vector<double> point(model.number_variables);
         Vector<double> jacobian(model.number_jacobian_nonzeros());
         for ([[maybe_unused]] size_t point_index: Range(number_random_points)) {
            for (size_t variable_index: Range(model.number_variables)) {
               const double value = base_point[variable_index];
               point[variable_index] = is_fixed[variable_index] ? value : std::clamp(value + 0.1 * std::max(1., std::abs(value)) *
                  distribution(generator), lower_bounds[variable_index], upper_bounds[variable_index]);
            }
            try {
               model.evaluate_jacobian(point, jacobian.data());
            }
            catch (const std::exception& exception) {
               WARNING << "Preprocessing: the linear constraints could not be detected (" << exception.what() << ")\n";
               return {};
            }
            // the derivatives w.r.t. the fixed variables are irrelevant: linearity is assessed in the reduced space
            for (size_t nonzero_index: Range(model.number_jacobian_nonzeros())) {
               if (!is_fixed[static_cast<size_t>(column_indices[nonzero_index])] &&
                     !are_close(jacobian[nonzero_index], base_jacobian[nonzero_index])) {
                  is_constant[static_cast<size_t>(row_indices[nonzero_index])] = false;
               }
            }
         }
         std::vector<size_t> linear_constraints;
         for (size_t constraint_index: model.get_nonlinear_constraints()) {
            if (is_constant[constraint_index]) {
               linear_constraints.push_back(constraint_index);
            }
         }
         return linear_constraints;
      }

      LinearRows linear_rows(const Model& model, const std::vector<bool>& is_linear, const Vector<double>& jacobian) {
         LinearRows rows(model.number_constraints);
         const auto row_indices = model.get_jacobian_row_indices();
         const auto column_indices = model.get_jacobian_column_indices();
         for (size_t nonzero_index: Range(model.number_jacobian_nonzeros())) {
            const size_t constraint_index = static_cast<size_t>(row_indices[nonzero_index]);
            if (is_linear[constraint_index]) {
               rows[constraint_index].push_back({static_cast<size_t>(column_indices[nonzero_index]), jacobian[nonzero_index]});
            }
         }
         return rows;
      }

      bool is_satisfied(double value, double lower_bound, double upper_bound) {
         const double tolerance = 1e-12 * std::max(1., std::abs(value));
         return lower_bound - tolerance <= value && value <= upper_bound + tolerance;
      }

      // turn the linear constraints with a single non-fixed variable into bounds, remove the satisfied linear constraints
      // without non-fixed variables, and fix the variables whose bounds become equal. Repeat until no variable gets fixed
      void eliminate_bound_constraints(const Model& model, const LinearRows& rows, const Vector<double>& base_point,
            const Vector<double>& constraints, std::vector<bool>& is_fixed, std::vector<bool>& is_removed,
            ModelReduction& reduction) {
         auto& lower_bounds = reduction.variables_lower_bounds;
         auto& upper_bounds = reduction.variables_upper_bounds;
         const auto& constraints_lower_bounds = model.get_constraints_lower_bounds();
         const auto& constraints_upper_bounds = model.get_constraints_upper_bounds();
         for (size_t round = 0; ; ++round) {
            for (size_t constraint_index: Range(model.number_constraints)) {
               const auto& terms = rows[constraint_index];
               if (is_removed[constraint_index] || terms.empty()) {
                  continue;
               }
               // the unique non-fixed variable x_i of the constraint
               std::optional<size_t> variable;
               bool several_variables = false;
               for (const LinearTerm& term: terms) {
                  if (!is_fixed[term.variable_index]) {
                     several_variables = several_variables || (variable.has_value() && *variable != term.variable_index);
                     variable = term.variable_index;
                  }
               }
               if (several_variables) {
                  continue;
               }
               // c_j(x) = a x_i + b, where the fixed variables are at their values (their lower bounds)
               double coefficient = 0.;
               double constant = constraints[constraint_index];
               for (const LinearTerm& term: terms) {
                  if (term.variable_index == variable) {
                     coefficient += term.coefficient;
                  }
                  else {
                     constant += term.coefficient * (lower_bounds[term.variable_index] - base_point[term.variable_index]);
                  }
               }
               if (variable.has_value()) {
                  constant -= coefficient * base_point[*variable];
               }
               if (!std::isfinite(coefficient) || !std::isfinite(constant)) {
                  continue;
               }
               // constant constraint: removed if satisfied, otherwise left to the solver, which will report the infeasibility
               if (coefficient == 0.) {
                  if (is_satisfied(constant, constraints_lower_bounds[constraint_index], constraints_upper_bounds[constraint_index])) {
                     is_removed[constraint_index] = true;
                     reduction.redundant_constraints.push_back(constraint_index);
                  }
                  continue;
               }
               const size_t variable_index = *variable;
               // cl <= a x_i + b <= cu  <=>  (cl - b)/a <= x_i <= (cu - b)/a if a > 0 (swapped if a < 0)
               double induced_lower_bound = (constraints_lower_bounds[constraint_index] - constant) / coefficient;
               double induced_upper_bound = (constraints_upper_bounds[constraint_index] - constant) / coefficient;
               if (coefficient < 0.) {
                  std::swap(induced_lower_bound, induced_upper_bound);
               }
               // inconsistent bounds: leave the constraint to the solver, which will report the infeasibility
               if (std::min(induced_upper_bound, upper_bounds[variable_index]) <
                     std::max(induced_lower_bound, lower_bounds[variable_index])) {
                  continue;
               }
               const size_t bound_constraint_index = reduction.bound_constraints.size();
               reduction.bound_constraints.push_back({constraint_index, variable_index, coefficient});
               is_removed[constraint_index] = true;
               if (lower_bounds[variable_index] < induced_lower_bound) {
                  lower_bounds[variable_index] = induced_lower_bound;
                  reduction.lower_bound_origins[variable_index] = bound_constraint_index;
               }
               if (induced_upper_bound < upper_bounds[variable_index]) {
                  upper_bounds[variable_index] = induced_upper_bound;
                  reduction.upper_bound_origins[variable_index] = bound_constraint_index;
               }
            }
            // the variables fixed in this round are eliminated in the next one
            bool new_fixed_variables = false;
            for (size_t variable_index: Range(model.number_variables)) {
               if (!is_fixed[variable_index] && lower_bounds[variable_index] == upper_bounds[variable_index]) {
                  is_fixed[variable_index] = new_fixed_variables = true;
                  reduction.fixed_variables.push_back(variable_index);
                  reduction.elimination_rounds.push_back(round + 1);
               }
            }
            if (!new_fixed_variables) {
               return;
            }
         }
      }
   } // namespace

   ModelReduction ModelReduction::detect(const Model& model, const Options& options) {
      ModelReduction reduction;
      reduction.variables_lower_bounds = model.get_variables_lower_bounds();
      reduction.variables_upper_bounds = model.get_variables_upper_bounds();
      reduction.lower_bound_origins.resize(model.number_variables);
      reduction.upper_bound_origins.resize(model.number_variables);

      // variables fixed in the model
      std::vector<bool> is_fixed(model.number_variables, false);
      for (size_t variable_index: Range(model.number_variables)) {
         if (is_finite(reduction.variables_lower_bounds[variable_index]) &&
               reduction.variables_lower_bounds[variable_index] == reduction.variables_upper_bounds[variable_index]) {
            is_fixed[variable_index] = true;
            reduction.fixed_variables.push_back(variable_index);
            reduction.elimination_rounds.push_back(0);
         }
      }

      std::vector<bool> is_removed(model.number_constraints, false);
      const bool bound_detection = options.get_bool("detect_bound_constraints");
      const bool linearity_detection = options.get_bool("detect_linear_constraints") &&
         !model.get_nonlinear_constraints().empty();
      if (0 < model.number_constraints && (bound_detection || linearity_detection)) {
         // base point: the initial point with the fixed variables at their values
         Vector<double> base_point(model.number_variables);
         model.initial_primal_point(base_point);
         for (size_t variable_index: reduction.fixed_variables) {
            base_point[variable_index] = reduction.variables_lower_bounds[variable_index];
         }
         Vector<double> constraints(model.number_constraints);
         Vector<double> jacobian(model.number_jacobian_nonzeros());
         bool evaluated = true;
         try {
            model.evaluate_constraints(base_point, constraints);
            model.evaluate_jacobian(base_point, jacobian.data());
         }
         catch (const std::exception& exception) {
            WARNING << "Preprocessing: the constraints could not be evaluated at the initial point (" << exception.what() << ")\n";
            evaluated = false;
         }
         if (evaluated && linearity_detection) {
            reduction.detected_linear_constraints = detect_linear_constraints(model, base_point, jacobian, is_fixed);
         }
         if (evaluated && bound_detection) {
            const LinearRows rows = linear_rows(model, reduction.linear_constraints(model), jacobian);
            eliminate_bound_constraints(model, rows, base_point, constraints, is_fixed, is_removed, reduction);
         }
      }

      for (size_t variable_index: Range(model.number_variables)) {
         if (!is_fixed[variable_index]) {
            reduction.free_variables.push_back(variable_index);
         }
      }
      for (size_t constraint_index: Range(model.number_constraints)) {
         if (!is_removed[constraint_index]) {
            reduction.kept_constraints.push_back(constraint_index);
         }
      }
      return reduction;
   }

   std::vector<bool> ModelReduction::linear_constraints(const Model& model) const {
      std::vector<bool> is_linear(model.number_constraints, true);
      for (size_t constraint_index: model.get_nonlinear_constraints()) {
         is_linear[constraint_index] = false;
      }
      for (size_t constraint_index: this->detected_linear_constraints) {
         is_linear[constraint_index] = true;
      }
      return is_linear;
   }

   bool ModelReduction::is_trivial() const {
      return this->fixed_variables.empty() && this->bound_constraints.empty() && this->redundant_constraints.empty() &&
         this->detected_linear_constraints.empty();
   }
} // namespace
