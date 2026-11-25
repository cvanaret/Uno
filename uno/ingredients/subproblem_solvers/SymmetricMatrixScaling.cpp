// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cmath>
#include <limits>
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
      else if (method == "matching") {
         return LinearSystemScalingMethod::MATCHING;
      }
      throw std::invalid_argument("The linear system scaling method " + method + " is unknown (none|ruiz|matching)");
   }

   void SymmetricMatrixScaling::initialize(size_t dimension, const uno_int* row_indices, const uno_int* column_indices,
         size_t number_nonzeros, int indexing) {
      this->dimension = dimension;
      this->previous_column_duals.clear(); // new pattern: no warm start
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
      else if (method == LinearSystemScalingMethod::MATCHING) {
         this->compute_matching_scaling(scaling);
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

   // symmetric Ruiz equilibration in the infinity norm, as in MC77 (JOB = 0, symmetric case): an initial pass, then at most
   // ruiz_max_iterations passes s_i <- s_i / sqrt(||row_i(S A S)||_inf). A pair of rows whose largest entries are the same
   // entry (i, j) is frozen: that entry is 1 after the pass, and so are both row norms from then on. Zero rows are left
   // unscaled. After each iteration (not the initial pass), stop if the row norms of the active rows before the pass are
   // all strictly within 1 +/- ruiz_tolerance
   void SymmetricMatrixScaling::compute_ruiz_scaling(std::vector<double>& scaling) const {
      std::vector<double> row_norms(this->dimension);
      std::vector<size_t> argmax(this->dimension); // column of the largest entry of each row
      std::vector<bool> active(this->dimension, true);
      for (size_t iteration = 0; iteration <= ruiz_max_iterations; ++iteration) {
         std::fill(row_norms.begin(), row_norms.end(), 0.);
         const auto update_row = [&](size_t row, size_t column, double value) {
            if (active[row] && row_norms[row] < value) {
               row_norms[row] = value;
               argmax[row] = column;
            }
         };
         for (size_t entry_index = 0; entry_index < this->entry_values.size(); ++entry_index) {
            const size_t row = this->entry_rows[entry_index], column = this->entry_columns[entry_index];
            const double value = std::abs(scaling[row] * this->entry_values[entry_index] * scaling[column]);
            update_row(column, row, value);
            update_row(row, column, value);
         }
         for (size_t index = 0; index < this->dimension; ++index) {
            if (active[index] && row_norms[index] == 0.) {
               active[index] = false; // zero row
            }
            else if (active[index]) {
               scaling[index] /= std::sqrt(row_norms[index]);
            }
         }
         double error = 0.;
         for (size_t index = 0; index < this->dimension; ++index) {
            if (active[index] && active[argmax[index]] && argmax[argmax[index]] == index) {
               active[index] = active[argmax[index]] = false;
            }
            else if (active[index]) {
               error = std::max(error, std::abs(1. - row_norms[index]));
            }
         }
         if (0 < iteration && error < ruiz_tolerance) {
            return;
         }
      }
   }

   namespace {
      constexpr double infinity = std::numeric_limits<double>::infinity();
   } // namespace

   // Symmetric scaling from a maximum product matching (Duff & Koster 2001, MC64 job 5; Duff & Pralet 2005).
   // 1. Bipartite graph of the full symmetric matrix A (both triangles, explicit zeros dropped): row i -- column j if
   //    a_ij != 0, with cost c_ij = log(max_k |a_kj|) - log|a_ij| >= 0. An assignment of minimum total cost maximizes the
   //    product of the absolute values of the matched entries.
   // 2. Minimum cost assignment with duals u (rows) and v (columns), see MinimumCostAssignment: u_i + v_j <= c_ij, with
   //    equality on the matching. The column duals of the previous matrix (same pattern) warm-start the next one, since
   //    consecutive KKT matrices are similar.
   // 3. Row and column scalings r_i = exp(u_i), c_j = exp(v_j) / max_k |a_kj|: |r_i a_ij c_j| <= 1, = 1 on the matching.
   //    Symmetric scaling s_i = sqrt(r_i c_i). Since A is symmetric, the transposed matching is also optimal, and
   //    complementary slackness makes the duals tight on both (i, j) and (j, i): |s_i a_ij s_j| =
   //    sqrt(|r_i a_ij c_j| |r_j a_ji c_i|) is at most 1, and 1 on the matched entries.
   // 4. Structurally singular matrices (Duff & Pralet 2005): the principal submatrix of the matched rows is structurally
   //    nonsingular, and its assignment (with the costs of A) yields the scaling of these indices; every other index is
   //    scaled such that the largest entry of its row, with the scaling of the matched indices, is 1.
   // These properties hold, up to rounding errors and the bounds [min_scaling, max_scaling], for the values the scaling is
   // computed from. Cost: one Dijkstra search per column left free by the cheap initial phases; on structurally singular
   // matrices, the searches from the columns that cannot be matched explore their whole component, and the submatrix is
   // solved without warm start.
   void SymmetricMatrixScaling::compute_matching_scaling(std::vector<double>& scaling) {
      const size_t dimension = this->dimension;
      std::vector<double> absolute_values; // |a_ij| of the edges
      std::vector<double> log_column_maximum; // log max_k |a_kj|
      BipartiteGraph graph = this->build_matching_graph(absolute_values, log_column_maximum);

      this->assignment.solve(graph, this->previous_column_duals);
      this->store_column_duals_for_warm_start(graph);
      const std::vector<size_t>& row_matches = this->assignment.get_row_matches();

      std::vector<bool> scaled_by_matching(dimension, false);
      bool structurally_singular = false;
      for (size_t index = 0; index < dimension; ++index) {
         if (row_matches[index] == MinimumCostAssignment::unmatched && !graph.is_empty(index)) {
            structurally_singular = true;
         }
      }
      if (!structurally_singular) {
         SymmetricMatrixScaling::scale_from_duals(this->assignment, log_column_maximum, {}, scaling, scaled_by_matching);
      }
      else {
         // principal submatrix of the matched rows, renumbered
         std::vector<size_t> submatrix_index(dimension, MinimumCostAssignment::unmatched), original_index;
         for (size_t index = 0; index < dimension; ++index) {
            if (row_matches[index] != MinimumCostAssignment::unmatched) {
               submatrix_index[index] = original_index.size();
               original_index.push_back(index);
            }
         }
         BipartiteGraph submatrix_graph;
         submatrix_graph.dimension = original_index.size();
         submatrix_graph.column_start.assign(1, 0);
         for (const size_t column: original_index) {
            for (size_t edge = graph.first_edge(column); edge < graph.end_edge(column); ++edge) {
               const size_t row = submatrix_index[graph.edge_rows[edge]];
               if (row != MinimumCostAssignment::unmatched) {
                  submatrix_graph.edge_rows.push_back(row);
                  submatrix_graph.edge_costs.push_back(graph.edge_costs[edge]);
               }
            }
            submatrix_graph.column_start.push_back(submatrix_graph.edge_rows.size());
         }
         MinimumCostAssignment submatrix_assignment;
         submatrix_assignment.solve(submatrix_graph);
         SymmetricMatrixScaling::scale_from_duals(submatrix_assignment, log_column_maximum, original_index, scaling,
            scaled_by_matching);
      }
      // the other indices: the largest entry of their row, with the scaling of the indices scaled by the matching, is 1
      for (size_t index = 0; index < dimension; ++index) {
         if (!scaled_by_matching[index]) {
            double row_maximum = 0.;
            for (size_t edge = graph.first_edge(index); edge < graph.end_edge(index); ++edge) {
               const size_t other_index = graph.edge_rows[edge];
               if (scaled_by_matching[other_index]) {
                  row_maximum = std::max(row_maximum, absolute_values[edge] * scaling[other_index]);
               }
            }
            scaling[index] = (0. < row_maximum) ? 1. / row_maximum : 1.;
         }
      }
   }

   // The column duals only increase during an assignment: they are shifted (the row duals are recomputed from them, so a
   // uniform shift does not change the warm start) such that their maximum over the nonempty columns is 0, which keeps
   // them small along a sequence of matrices (and the rounding errors in the scaling with them)
   void SymmetricMatrixScaling::store_column_duals_for_warm_start(const BipartiteGraph& graph) {
      this->previous_column_duals = this->assignment.get_column_duals();
      double maximum = -infinity;
      for (size_t column = 0; column < graph.dimension; ++column) {
         if (!graph.is_empty(column)) {
            maximum = std::max(maximum, this->previous_column_duals[column]);
         }
      }
      if (std::isfinite(maximum)) {
         for (double& dual: this->previous_column_duals) {
            dual -= maximum;
         }
      }
   }

   // bipartite graph of the full symmetric pattern (column j: the rows i with a_ij != 0) with the costs
   // c_ij = log(max_k |a_kj|) - log|a_ij|
   BipartiteGraph SymmetricMatrixScaling::build_matching_graph(std::vector<double>& absolute_values,
         std::vector<double>& log_column_maximum) const {
      const size_t dimension = this->dimension;
      BipartiteGraph graph;
      graph.dimension = dimension;
      // number of edges of each column, then their positions
      graph.column_start.assign(dimension + 1, 0);
      for (size_t entry_index = 0; entry_index < this->entry_values.size(); ++entry_index) {
         if (this->entry_values[entry_index] != 0.) {
            ++graph.column_start[this->entry_columns[entry_index] + 1];
            if (this->entry_rows[entry_index] != this->entry_columns[entry_index]) {
               ++graph.column_start[this->entry_rows[entry_index] + 1];
            }
         }
      }
      for (size_t column = 0; column < dimension; ++column) {
         graph.column_start[column + 1] += graph.column_start[column];
      }
      const size_t number_edges = graph.column_start[dimension];
      graph.edge_rows.resize(number_edges);
      graph.edge_costs.resize(number_edges);
      absolute_values.resize(number_edges);
      std::vector<size_t> next_position(graph.column_start.begin(), graph.column_start.end() - 1);
      const auto add_edge = [&](size_t row, size_t column, double absolute_value) {
         const size_t edge = next_position[column]++;
         graph.edge_rows[edge] = row;
         absolute_values[edge] = absolute_value;
      };
      for (size_t entry_index = 0; entry_index < this->entry_values.size(); ++entry_index) {
         const double absolute_value = std::abs(this->entry_values[entry_index]);
         if (absolute_value != 0.) {
            const size_t row = this->entry_rows[entry_index], column = this->entry_columns[entry_index];
            add_edge(row, column, absolute_value);
            if (row != column) {
               add_edge(column, row, absolute_value);
            }
         }
      }
      // costs
      log_column_maximum.assign(dimension, -infinity);
      for (size_t column = 0; column < dimension; ++column) {
         for (size_t edge = graph.first_edge(column); edge < graph.end_edge(column); ++edge) {
            graph.edge_costs[edge] = std::log(absolute_values[edge]);
            log_column_maximum[column] = std::max(log_column_maximum[column], graph.edge_costs[edge]);
         }
         for (size_t edge = graph.first_edge(column); edge < graph.end_edge(column); ++edge) {
            graph.edge_costs[edge] = log_column_maximum[column] - graph.edge_costs[edge];
         }
      }
      return graph;
   }

   // s_i = sqrt(r_i c_i) = exp((u_i + v_i - log max_k |a_ki|) / 2) of the indices matched as a row and as a column. The
   // assignment may be that of a principal submatrix, whose indices are mapped by original_index (empty: identity)
   void SymmetricMatrixScaling::scale_from_duals(const MinimumCostAssignment& assignment,
         const std::vector<double>& log_column_maximum, const std::vector<size_t>& original_index,
         std::vector<double>& scaling, std::vector<bool>& scaled_by_matching) {
      const std::vector<double>& row_duals = assignment.get_row_duals();
      const std::vector<double>& column_duals = assignment.get_column_duals();
      const std::vector<size_t>& row_matches = assignment.get_row_matches();
      const std::vector<size_t>& column_matches = assignment.get_column_matches();
      for (size_t index = 0; index < row_duals.size(); ++index) {
         if (row_matches[index] != MinimumCostAssignment::unmatched && column_matches[index] != MinimumCostAssignment::unmatched) {
            const size_t matrix_index = original_index.empty() ? index : original_index[index];
            scaling[matrix_index] = std::exp(0.5 * (row_duals[index] + column_duals[index] - log_column_maximum[matrix_index]));
            scaled_by_matching[matrix_index] = true;
         }
      }
   }
} // namespace
