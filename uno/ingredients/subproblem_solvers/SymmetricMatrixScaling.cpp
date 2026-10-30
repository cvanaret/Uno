// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <utility>
#include "SymmetricMatrixScaling.hpp"

namespace uno {
   LinearSystemScalingMethod linear_system_scaling_method_from_string(const std::string& method) {
      if (method == "none") {
         return LinearSystemScalingMethod::NONE;
      }
      else if (method == "ruiz") {
         return LinearSystemScalingMethod::RUIZ;
      }
      throw std::invalid_argument("The linear system scaling method " + method + " is unknown (none|ruiz)");
   }

   void SymmetricMatrixScaling::initialize(size_t dimension, const uno_int* row_indices, const uno_int* column_indices,
         size_t number_nonzeros, int indexing) {
      this->dimension = dimension;
      this->coo_rows.resize(number_nonzeros);
      this->coo_columns.resize(number_nonzeros);
      this->coo_to_entry.resize(number_nonzeros);
      // lower triangle, 0-based
      for (size_t nonzero_index = 0; nonzero_index < number_nonzeros; ++nonzero_index) {
         const auto row = static_cast<size_t>(row_indices[nonzero_index] - indexing);
         const auto column = static_cast<size_t>(column_indices[nonzero_index] - indexing);
         if (dimension <= row || dimension <= column) {
            throw std::out_of_range("SymmetricMatrixScaling: index out of range");
         }
         this->coo_rows[nonzero_index] = row;
         this->coo_columns[nonzero_index] = column;
      }
      // compress the duplicates: sort the COO entries by (column, row) of the lower triangle
      std::vector<size_t> order(number_nonzeros);
      std::iota(order.begin(), order.end(), size_t(0));
      const auto lower = [&](size_t index) {
         const size_t row = this->coo_rows[index], column = this->coo_columns[index];
         return std::make_pair(std::min(row, column), std::max(row, column)); // (column, row) of the lower triangle
      };
      std::sort(order.begin(), order.end(), [&](size_t index1, size_t index2) { return lower(index1) < lower(index2); });
      this->entry_rows.clear();
      this->entry_columns.clear();
      for (size_t position = 0; position < number_nonzeros; ++position) {
         const auto [column, row] = lower(order[position]);
         if (position == 0 || lower(order[position - 1]) != std::make_pair(column, row)) {
            this->entry_rows.push_back(row);
            this->entry_columns.push_back(column);
         }
         this->coo_to_entry[order[position]] = this->entry_rows.size() - 1;
      }
      const size_t number_entries = this->entry_rows.size();
      this->entry_values.resize(number_entries);

      // full (both triangles) column-wise adjacency
      std::vector<size_t> column_count(dimension, 0);
      for (size_t entry_index = 0; entry_index < number_entries; ++entry_index) {
         ++column_count[this->entry_columns[entry_index]];
         if (this->entry_rows[entry_index] != this->entry_columns[entry_index]) {
            ++column_count[this->entry_rows[entry_index]];
         }
      }
      this->initialized = true;
   }

   void SymmetricMatrixScaling::compute(LinearSystemScalingMethod method, const double* coo_values, std::vector<double>& scaling) {
      if (!this->initialized) {
         throw std::logic_error("SymmetricMatrixScaling::compute: the pattern was not initialized");
      }
      scaling.assign(this->dimension, 1.);
      this->assemble_entry_values(coo_values);
      // non-finite values: the scaling would be meaningless, leave the matrix unscaled
      if (!std::all_of(this->entry_values.begin(), this->entry_values.end(), [](double value) { return std::isfinite(value); })) {
         return;
      }
      if (method == LinearSystemScalingMethod::RUIZ) {
         this->compute_ruiz_scaling(scaling);
      }
      for (double& factor: scaling) {
         factor = std::clamp(factor, min_scaling, max_scaling);
      }
   }

   void SymmetricMatrixScaling::compute_entry_factors(const std::vector<double>& scaling, std::vector<double>& entry_factors) const {
      const size_t number_nonzeros = this->coo_rows.size();
      entry_factors.resize(number_nonzeros);
      for (size_t nonzero_index = 0; nonzero_index < number_nonzeros; ++nonzero_index) {
         entry_factors[nonzero_index] = scaling[this->coo_rows[nonzero_index]] * scaling[this->coo_columns[nonzero_index]];
      }
   }

   void SymmetricMatrixScaling::assemble_entry_values(const double* coo_values) {
      std::fill(this->entry_values.begin(), this->entry_values.end(), 0.);
      for (size_t nonzero_index = 0; nonzero_index < this->coo_to_entry.size(); ++nonzero_index) {
         this->entry_values[this->coo_to_entry[nonzero_index]] += coo_values[nonzero_index];
      }
   }

   // symmetric Ruiz equilibration: s_i <- s_i / sqrt(||row_i(S A S)||_inf) until all the row norms are close to 1
   void SymmetricMatrixScaling::compute_ruiz_scaling(std::vector<double>& scaling) const {
      std::vector<double> row_norms(this->dimension);
      for (size_t iteration = 0; iteration < ruiz_max_iterations; ++iteration) {
         std::fill(row_norms.begin(), row_norms.end(), 0.);
         for (size_t entry_index = 0; entry_index < this->entry_values.size(); ++entry_index) {
            const size_t row = this->entry_rows[entry_index], column = this->entry_columns[entry_index];
            const double value = std::abs(scaling[row] * this->entry_values[entry_index] * scaling[column]);
            row_norms[row] = std::max(row_norms[row], value);
            row_norms[column] = std::max(row_norms[column], value);
         }
         bool converged = true;
         for (size_t index = 0; index < this->dimension; ++index) {
            if (0. < row_norms[index]) {
               converged = converged && std::abs(1. - row_norms[index]) <= ruiz_tolerance;
               scaling[index] /= std::sqrt(row_norms[index]);
            }
         }
         if (converged) {
            return;
         }
      }
   }
} // namespace
