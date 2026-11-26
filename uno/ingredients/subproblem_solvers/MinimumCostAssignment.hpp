// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MINIMUMCOSTASSIGNMENT_H
#define UNO_MINIMUMCOSTASSIGNMENT_H

#include <cstddef>
#include <limits>
#include <vector>

namespace uno {
   class IndexedMinHeap;

   // sparse bipartite graph with n rows and n columns, in compressed column format: the edges of column j are
   // column_start[j], ..., column_start[j + 1] - 1, edge e joins row edge_rows[e] and column j with cost edge_costs[e]
   struct BipartiteGraph {
      size_t dimension{0};
      std::vector<size_t> column_start{0}; // size dimension + 1
      std::vector<size_t> edge_rows{};
      std::vector<double> edge_costs{};

      [[nodiscard]] size_t first_edge(size_t column) const { return this->column_start[column]; }
      [[nodiscard]] size_t end_edge(size_t column) const { return this->column_start[column + 1]; }
      [[nodiscard]] bool is_empty(size_t column) const { return this->first_edge(column) == this->end_edge(column); }
   };

   // Minimum cost assignment (linear sum assignment) on a sparse bipartite graph with nonnegative costs: among the
   // matchings of maximum cardinality, find one with minimum total cost. Sparse Hungarian algorithm (primal-dual):
   // the dual variables (u_i of the rows, v_j of the columns) remain feasible, i.e. the reduced costs
   // c_ij - u_i - v_j are nonnegative, and the matched edges are tight (zero reduced cost). When the matching has maximum
   // cardinality, complementary slackness certifies its optimality.
   // 1. initial duals: column duals given by the caller (warm start, e.g. from a similar previous graph) or 0, then
   //    u_i = min_j (c_ij - v_j) and v_j = min_i (c_ij - u_i);
   // 2. cheap initial matching on the tight edges: greedy, then augmenting paths of length 2;
   // 3. augmenting row reduction (Jonker & Volgenant 1987): free columns bid for their cheapest row;
   // 4. one shortest augmenting path (Dijkstra on the reduced costs) from each remaining free column. The rows are
   //    labeled with their distance from the free column; the rows at the current minimum distance are processed from a
   //    stack, the others wait in an indexed heap. The search stops as soon as no labeled row is closer than the best
   //    free row reached so far, and the labels that cannot improve on it are discarded.
   // References: Duff & Koster 2001 (MC64), Jonker & Volgenant 1987, Burkard, Dell'Amico & Martello 2012 (Assignment
   // problems, chapter 4)
   class MinimumCostAssignment {
   public:
      static constexpr size_t unmatched{std::numeric_limits<size_t>::max()};

      MinimumCostAssignment() = default;

      // compute an optimal assignment and duals. initial_column_duals: warm start, ignored if its size is not
      // graph.dimension (non-finite entries are replaced by 0)
      void solve(const BipartiteGraph& graph, const std::vector<double>& initial_column_duals = {});

      // results of the last solve
      [[nodiscard]] const std::vector<double>& get_row_duals() const { return this->row_duals; }
      [[nodiscard]] const std::vector<double>& get_column_duals() const { return this->column_duals; }
      // column matched to each row, or unmatched
      [[nodiscard]] const std::vector<size_t>& get_row_matches() const { return this->row_match; }
      // row matched to each column, or unmatched
      [[nodiscard]] const std::vector<size_t>& get_column_matches() const { return this->column_match; }

   protected:
      std::vector<double> row_duals{};
      std::vector<double> column_duals{};
      std::vector<size_t> row_match{};
      std::vector<size_t> column_match{};

      // workspace of the shortest path searches
      std::vector<double> row_distance{};
      std::vector<double> column_distance{};
      std::vector<size_t> row_predecessor{}; // column from which the row was labeled
      std::vector<char> row_finalized{};
      std::vector<size_t> labeled_rows{};
      std::vector<size_t> finalized_rows{};
      std::vector<size_t> finalized_columns{};
      std::vector<size_t> rows_at_current_distance{};

      [[nodiscard]] double reduced_cost(const BipartiteGraph& graph, size_t edge, size_t column) const;
      void match(size_t row, size_t column);
      void initialize_duals(const BipartiteGraph& graph, const std::vector<double>& initial_column_duals);
      void match_greedily(const BipartiteGraph& graph);
      void augment_along_paths_of_length_two(const BipartiteGraph& graph);
      void reduce_augmenting_rows(const BipartiteGraph& graph);
      void augment_along_shortest_paths(const BipartiteGraph& graph);
      [[nodiscard]] size_t find_shortest_augmenting_path(const BipartiteGraph& graph, size_t free_column, IndexedMinHeap& heap,
         double& path_length);
   };
} // namespace

#endif // UNO_MINIMUMCOSTASSIGNMENT_H
