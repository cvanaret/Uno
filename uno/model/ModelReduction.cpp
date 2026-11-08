// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cmath>
#include <exception>
#include "ModelReduction.hpp"
#include "Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/Range.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   namespace {
      // constraints whose Jacobian row has no entry in a non-fixed column
      std::vector<size_t> constraints_of_fixed_variables(const Model& model, const std::vector<bool>& is_fixed) {
         std::vector<bool> depends_on_free_variables(model.number_constraints, false);
         const auto row_indices = model.get_jacobian_row_indices();
         const auto column_indices = model.get_jacobian_column_indices();
         for (size_t nonzero_index: Range(model.number_jacobian_nonzeros())) {
            if (!is_fixed[static_cast<size_t>(column_indices[nonzero_index])]) {
               depends_on_free_variables[static_cast<size_t>(row_indices[nonzero_index])] = true;
            }
         }
         std::vector<size_t> constraints;
         for (size_t constraint_index: Range(model.number_constraints)) {
            if (!depends_on_free_variables[constraint_index]) {
               constraints.push_back(constraint_index);
            }
         }
         return constraints;
      }

      // satisfied up to roundoff, but never violated by more than feasibility_tolerance
      bool is_satisfied(double value, double lower_bound, double upper_bound, double feasibility_tolerance) {
         const double tolerance = std::min(1e-12 * std::max(1., std::abs(value)), feasibility_tolerance);
         return std::isfinite(value) && lower_bound - tolerance <= value && value <= upper_bound + tolerance;
      }
   } // namespace

   ModelReduction ModelReduction::detect(const Model& model, double feasibility_tolerance) {
      ModelReduction reduction;
      const auto& lower_bounds = model.get_variables_lower_bounds();
      const auto& upper_bounds = model.get_variables_upper_bounds();
      std::vector<bool> is_fixed(model.number_variables);
      for (size_t variable_index: Range(model.number_variables)) {
         is_fixed[variable_index] = is_finite(lower_bounds[variable_index]) &&
            lower_bounds[variable_index] == upper_bounds[variable_index];
         (is_fixed[variable_index] ? reduction.fixed_variables : reduction.free_variables).push_back(variable_index);
      }

      // the constraints of fixed variables only are evaluated at the initial point, with the fixed variables at their values
      std::vector<bool> is_removed(model.number_constraints, false);
      const std::vector<size_t> constant_constraints = reduction.fixed_variables.empty() ? std::vector<size_t>{} :
         constraints_of_fixed_variables(model, is_fixed);
      if (!constant_constraints.empty()) {
         Vector<double> point(model.number_variables);
         model.initial_primal_point(point);
         for (size_t variable_index: reduction.fixed_variables) {
            point[variable_index] = lower_bounds[variable_index];
         }
         Vector<double> constraints(model.number_constraints);
         try {
            model.evaluate_constraints(point, constraints);
            const auto& constraints_lower_bounds = model.get_constraints_lower_bounds();
            const auto& constraints_upper_bounds = model.get_constraints_upper_bounds();
            for (size_t constraint_index: constant_constraints) {
               is_removed[constraint_index] = is_satisfied(constraints[constraint_index], constraints_lower_bounds[constraint_index],
                  constraints_upper_bounds[constraint_index], feasibility_tolerance);
            }
         }
         catch (const std::exception& exception) {
            WARNING << "Preprocessing: the constraints of fixed variables only could not be evaluated (" << exception.what() << ")\n";
         }
      }
      for (size_t constraint_index: Range(model.number_constraints)) {
         (is_removed[constraint_index] ? reduction.removed_constraints : reduction.kept_constraints).push_back(constraint_index);
      }
      return reduction;
   }

   bool ModelReduction::is_trivial() const {
      return this->fixed_variables.empty();
   }
} // namespace
