// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <stdexcept>
#include <string>
#include "QuadraticProgram.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   // out-of-line anchor for the vtable
   QuadraticProgram::~QuadraticProgram() = default;

   namespace {
      void validate_coo_matrix(const std::string& name, const Vector<uno_int>& row_indices,
            const Vector<uno_int>& column_indices, const Vector<double>& values, size_t number_rows, size_t number_columns) {
         if (row_indices.size() != values.size() || column_indices.size() != values.size()) {
            throw std::invalid_argument("QuadraticProgram: the " + name + " COO arrays have inconsistent lengths (" +
               std::to_string(row_indices.size()) + ", " + std::to_string(column_indices.size()) + ", " +
               std::to_string(values.size()) + ")");
         }
         for (size_t nonzero_index: Range(values.size())) {
            const uno_int row_index = row_indices[nonzero_index];
            const uno_int column_index = column_indices[nonzero_index];
            if (row_index < 0 || number_rows <= static_cast<size_t>(row_index) ||
                  column_index < 0 || number_columns <= static_cast<size_t>(column_index)) {
               throw std::invalid_argument("QuadraticProgram: the " + name + " entry (" + std::to_string(row_index) + ", " +
                  std::to_string(column_index) + ") is out of range [0, " + std::to_string(number_rows) + ") x [0, " +
                  std::to_string(number_columns) + ")");
            }
            if (!std::isfinite(values[nonzero_index])) {
               throw std::invalid_argument("QuadraticProgram: the " + name + " has a non-finite entry");
            }
         }
      }

      void validate_bounds(const std::string& name, const std::vector<double>& lower_bounds,
            const std::vector<double>& upper_bounds, size_t dimension) {
         if (lower_bounds.size() != dimension || upper_bounds.size() != dimension) {
            throw std::invalid_argument("QuadraticProgram: the " + name + " bounds should have size " + std::to_string(dimension));
         }
         for (size_t index: Range(dimension)) {
            if (std::isnan(lower_bounds[index]) || std::isnan(upper_bounds[index])) {
               throw std::invalid_argument("QuadraticProgram: the " + name + " bounds contain NaN");
            }
         }
      }
   } // namespace

   void QuadraticProgram::validate_coo_data(const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values,
         const Vector<uno_int>& hessian_row_indices, const Vector<uno_int>& hessian_column_indices,
         const Vector<double>& hessian_values,
         const std::vector<double>& variables_lower_bounds, const std::vector<double>& variables_upper_bounds,
         const std::vector<double>& constraints_lower_bounds, const std::vector<double>& constraints_upper_bounds) {
      const size_t number_variables = linear_objective.size();
      const size_t number_constraints = constraints_lower_bounds.size();
      for (size_t variable_index: Range(number_variables)) {
         if (!std::isfinite(linear_objective[variable_index])) {
            throw std::invalid_argument("QuadraticProgram: the linear objective has a non-finite entry");
         }
      }
      validate_coo_matrix("Jacobian", jacobian_row_indices, jacobian_column_indices, jacobian_values, number_constraints,
         number_variables);
      validate_coo_matrix("Hessian", hessian_row_indices, hessian_column_indices, hessian_values, number_variables,
         number_variables);
      validate_bounds("variable", variables_lower_bounds, variables_upper_bounds, number_variables);
      validate_bounds("constraint", constraints_lower_bounds, constraints_upper_bounds, number_constraints);
   }
} // namespace
