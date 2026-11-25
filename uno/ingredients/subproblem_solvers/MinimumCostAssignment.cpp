// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cmath>
#include <utility>
#include "MinimumCostAssignment.hpp"
#include "tools/IndexedMinHeap.hpp"

namespace uno {
   namespace {
      constexpr double infinity = std::numeric_limits<double>::infinity();
   } // namespace

   void MinimumCostAssignment::solve(const BipartiteGraph& graph, const std::vector<double>& initial_column_duals) {
      this->initialize_duals(graph, initial_column_duals);
      this->match_greedily(graph);
      this->augment_along_paths_of_length_two(graph);
      this->reduce_augmenting_rows(graph);
      this->augment_along_shortest_paths(graph);
   }

   // the reduced cost is evaluated as (c_ij - u_i) - v_j, the order in which the duals are computed, so that the
   // minima of the initial duals are exactly 0 (the edges are tight without a tolerance)
   double MinimumCostAssignment::reduced_cost(const BipartiteGraph& graph, size_t edge, size_t column) const {
      return (graph.edge_costs[edge] - this->row_duals[graph.edge_rows[edge]]) - this->column_duals[column];
   }

   void MinimumCostAssignment::match(size_t row, size_t column) {
      this->row_match[row] = column;
      this->column_match[column] = row;
   }

   // u_i = min_j (c_ij - v_j), with the initial v (0 for a row without edges)
   void MinimumCostAssignment::initialize_duals(const BipartiteGraph& graph, const std::vector<double>& initial_column_duals) {
      const size_t dimension = graph.dimension;
      this->row_duals.assign(dimension, infinity);
      this->column_duals.assign(dimension, 0.);
      if (initial_column_duals.size() == dimension) {
         for (size_t column = 0; column < dimension; ++column) {
            if (std::isfinite(initial_column_duals[column])) {
               this->column_duals[column] = initial_column_duals[column];
            }
         }
      }
      for (size_t column = 0; column < dimension; ++column) {
         for (size_t edge = graph.first_edge(column); edge < graph.end_edge(column); ++edge) {
            double& row_dual = this->row_duals[graph.edge_rows[edge]];
            row_dual = std::min(row_dual, graph.edge_costs[edge] - this->column_duals[column]);
         }
      }
      for (double& row_dual: this->row_duals) {
         if (row_dual == infinity) {
            row_dual = 0.;
         }
      }
   }

   // v_j = min_i (c_ij - u_i), then each column takes the first free row of its tight edges
   void MinimumCostAssignment::match_greedily(const BipartiteGraph& graph) {
      const size_t dimension = graph.dimension;
      this->row_match.assign(dimension, unmatched);
      this->column_match.assign(dimension, unmatched);
      for (size_t column = 0; column < dimension; ++column) {
         if (graph.is_empty(column)) {
            continue;
         }
         double minimum = infinity;
         for (size_t edge = graph.first_edge(column); edge < graph.end_edge(column); ++edge) {
            minimum = std::min(minimum, graph.edge_costs[edge] - this->row_duals[graph.edge_rows[edge]]);
         }
         this->column_duals[column] = minimum;
         for (size_t edge = graph.first_edge(column); edge < graph.end_edge(column); ++edge) {
            const size_t row = graph.edge_rows[edge];
            if (this->row_match[row] == unmatched && this->reduced_cost(graph, edge, column) <= 0.) {
               this->match(row, column);
               break;
            }
         }
      }
   }

   // augmenting paths of length 2 on the tight edges: a free column j takes the row i of a tight edge if the column k
   // matched to i can move to a free row along a tight edge. The matched rows remain matched in this phase, so the search
   // for a free row in column k resumes where the previous search in column k stopped: each column is scanned once
   void MinimumCostAssignment::augment_along_paths_of_length_two(const BipartiteGraph& graph) {
      const size_t dimension = graph.dimension;
      std::vector<size_t> next_edge_to_scan(graph.column_start.begin(), graph.column_start.end() - 1);
      for (size_t column = 0; column < dimension; ++column) {
         for (size_t edge = graph.first_edge(column); edge < graph.end_edge(column) && this->column_match[column] == unmatched;
               ++edge) {
            if (0. < this->reduced_cost(graph, edge, column)) {
               continue;
            }
            const size_t row = graph.edge_rows[edge];
            const size_t other_column = this->row_match[row];
            if (other_column == unmatched) { // cannot happen after the greedy matching, kept for robustness
               this->match(row, column);
               break;
            }
            for (size_t& other_edge = next_edge_to_scan[other_column]; other_edge < graph.end_edge(other_column); ++other_edge) {
               const size_t other_row = graph.edge_rows[other_edge];
               if (this->row_match[other_row] == unmatched && this->reduced_cost(graph, other_edge, other_column) <= 0.) {
                  this->match(other_row, other_column);
                  this->match(row, column);
                  ++other_edge;
                  break;
               }
            }
         }
      }
   }

   // augmenting row reduction (Jonker & Volgenant 1987): a free column j takes its row i of smallest reduced cost r1,
   // and u_i is decreased by r2 - r1 (r2: second smallest reduced cost in column j), which makes (i, j) tight while
   // keeping the duals feasible; the other columns of row i become less attractive. The column previously matched to i
   // becomes free: it bids again in this pass if u_i strictly decreased, otherwise in the next pass. Two passes, each
   // with a bounded number of bids (the bids may cycle on ties)
   void MinimumCostAssignment::reduce_augmenting_rows(const BipartiteGraph& graph) {
      const size_t dimension = graph.dimension;
      std::vector<size_t> free_columns;
      for (size_t column = 0; column < dimension; ++column) {
         if (this->column_match[column] == unmatched && !graph.is_empty(column)) {
            free_columns.push_back(column);
         }
      }
      constexpr size_t number_passes = 2;
      for (size_t pass = 0; pass < number_passes && !free_columns.empty(); ++pass) {
         std::vector<size_t> next_free_columns;
         size_t position = 0;
         size_t remaining_bids = 4 * free_columns.size() + dimension;
         while (position < free_columns.size() && 0 < remaining_bids) {
            --remaining_bids;
            const size_t column = free_columns[position++];
            // the two edges of smallest reduced costs
            size_t best_edge = unmatched, second_best_edge = unmatched;
            double best_reduced_cost = infinity, second_best_reduced_cost = infinity;
            for (size_t edge = graph.first_edge(column); edge < graph.end_edge(column); ++edge) {
               const double reduced_cost = this->reduced_cost(graph, edge, column);
               if (reduced_cost < best_reduced_cost) {
                  second_best_reduced_cost = best_reduced_cost;
                  second_best_edge = best_edge;
                  best_reduced_cost = reduced_cost;
                  best_edge = edge;
               }
               else if (reduced_cost < second_best_reduced_cost) {
                  second_best_reduced_cost = reduced_cost;
                  second_best_edge = edge;
               }
            }
            // tie: prefer a free row (the augmentation then ends)
            size_t edge = best_edge;
            if (second_best_edge != unmatched && second_best_reduced_cost == best_reduced_cost &&
                  this->row_match[graph.edge_rows[best_edge]] != unmatched &&
                  this->row_match[graph.edge_rows[second_best_edge]] == unmatched) {
               edge = second_best_edge;
            }
            const size_t row = graph.edge_rows[edge];
            const double dual_decrease = (second_best_edge != unmatched && edge == best_edge) ?
               second_best_reduced_cost - best_reduced_cost : 0.;
            this->row_duals[row] -= dual_decrease;
            this->column_duals[column] = graph.edge_costs[edge] - this->row_duals[row]; // (row, column) becomes tight
            const size_t previous_column = this->row_match[row];
            this->match(row, column);
            if (previous_column != unmatched) {
               this->column_match[previous_column] = unmatched;
               if (0. < dual_decrease) {
                  free_columns.push_back(previous_column);
               }
               else {
                  next_free_columns.push_back(previous_column);
               }
            }
         }
         // the columns that did not bid in this pass
         next_free_columns.insert(next_free_columns.end(), free_columns.begin() + static_cast<std::ptrdiff_t>(position),
            free_columns.end());
         free_columns = std::move(next_free_columns);
      }
   }

   // one shortest augmenting path from each free column. After the search of a path of length D (distances measured with
   // the reduced costs), the duals of the finalized rows (distance d_i < D) and columns (distance d_j) are updated:
   // u_i += d_i - D, v_j += D - d_j. The reduced costs remain nonnegative and become 0 along the path, which is then
   // augmented. Without a path (no free row reachable), the column remains unmatched
   void MinimumCostAssignment::augment_along_shortest_paths(const BipartiteGraph& graph) {
      const size_t dimension = graph.dimension;
      this->row_distance.assign(dimension, infinity);
      this->column_distance.assign(dimension, infinity);
      this->row_predecessor.assign(dimension, unmatched);
      this->row_finalized.assign(dimension, 0);
      IndexedMinHeap heap(dimension, this->row_distance);
      for (size_t free_column = 0; free_column < dimension; ++free_column) {
         if (this->column_match[free_column] != unmatched || graph.is_empty(free_column)) {
            continue;
         }
         double path_length = infinity;
         const size_t free_row = this->find_shortest_augmenting_path(graph, free_column, heap, path_length);
         if (free_row != unmatched) {
            for (const size_t row: this->finalized_rows) {
               this->row_duals[row] += this->row_distance[row] - path_length;
            }
            for (const size_t column: this->finalized_columns) {
               this->column_duals[column] += path_length - this->column_distance[column];
            }
            // augmentation: flip the alternating path from the free row back to the free column
            size_t row = free_row;
            while (true) {
               const size_t column = this->row_predecessor[row];
               const size_t previous_row = this->column_match[column];
               this->match(row, column);
               if (column == free_column) {
                  break;
               }
               row = previous_row;
            }
         }
         // reset the labels of this search
         for (const size_t row: this->labeled_rows) {
            this->row_distance[row] = infinity;
            this->row_finalized[row] = 0;
            this->row_predecessor[row] = unmatched;
         }
         for (const size_t column: this->finalized_columns) {
            this->column_distance[column] = infinity;
         }
         this->labeled_rows.clear();
         this->finalized_rows.clear();
         this->finalized_columns.clear();
         this->rows_at_current_distance.clear();
         heap.clear();
      }
   }

   // Dijkstra from free_column on the reduced costs. Each finalized row leads to its matched column, whose edges label
   // further rows. Returns the free row at the end of a shortest augmenting path (and its length), or unmatched
   size_t MinimumCostAssignment::find_shortest_augmenting_path(const BipartiteGraph& graph, size_t free_column,
         IndexedMinHeap& heap, double& path_length) {
      size_t best_free_row = unmatched; // closest free row labeled so far, at distance path_length
      path_length = infinity;
      double current_distance = 0.;
      const auto scan_column = [&](size_t column, double distance) {
         this->column_distance[column] = distance;
         this->finalized_columns.push_back(column);
         for (size_t edge = graph.first_edge(column); edge < graph.end_edge(column); ++edge) {
            const size_t row = graph.edge_rows[edge];
            if (this->row_finalized[row]) {
               continue;
            }
            // the reduced costs are nonnegative up to rounding errors
            const double distance_via_column = distance + std::max(0., this->reduced_cost(graph, edge, column));
            // keep only the labels that improve on the row and on the best free row (in particular, best_free_row remains
            // the closest free row)
            if (path_length <= distance_via_column || this->row_distance[row] <= distance_via_column) {
               continue;
            }
            if (this->row_distance[row] == infinity) {
               this->labeled_rows.push_back(row);
            }
            this->row_distance[row] = distance_via_column;
            this->row_predecessor[row] = column;
            if (this->row_match[row] == unmatched) {
               best_free_row = row;
               path_length = distance_via_column;
            }
            else if (distance_via_column <= current_distance) {
               // at the current minimum distance: processed from the stack, without the heap
               if (heap.contains(row)) {
                  heap.remove(row);
               }
               this->rows_at_current_distance.push_back(row);
            }
            else {
               heap.push_or_decrease(row);
            }
         }
      };

      scan_column(free_column, 0.);
      while (true) {
         if (this->rows_at_current_distance.empty()) {
            if (heap.empty() || path_length <= this->row_distance[heap.top()]) {
               break;
            }
            // move all the rows at the new minimum distance to the stack
            current_distance = this->row_distance[heap.top()];
            while (!heap.empty() && this->row_distance[heap.top()] <= current_distance) {
               this->rows_at_current_distance.push_back(heap.pop());
            }
         }
         const size_t row = this->rows_at_current_distance.back();
         this->rows_at_current_distance.pop_back();
         if (this->row_finalized[row]) {
            continue; // stacked twice
         }
         if (path_length <= this->row_distance[row]) {
            break; // no labeled row is closer than the best free row
         }
         this->row_finalized[row] = 1;
         this->finalized_rows.push_back(row);
         scan_column(this->row_match[row], this->row_distance[row]);
      }
      return best_free_row;
   }
} // namespace
