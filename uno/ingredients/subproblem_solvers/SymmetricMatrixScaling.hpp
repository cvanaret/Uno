// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICMATRIXSCALING_H
#define UNO_SYMMETRICMATRIXSCALING_H

#include <cstddef>
#include <string>
#include <vector>
#include "MinimumCostAssignment.hpp"
#include "../interfaces/C/uno_int.h"

namespace uno {
   enum class LinearSystemScalingMethod {
      NONE,
      RUIZ, // symmetric Ruiz infinity-norm equilibration (Ruiz 2001; Knight, Ruiz & Uçar 2014)
      // symmetric scaling from a maximum product matching (Duff & Koster 2001, MC64 job 5; symmetrized as in Duff &
      // Pralet 2005, HSL_MC64): all the entries of S A S are at most 1 in absolute value, the matched entries are 1 (up
      // to rounding errors and the bounds on the scaling factors)
      MATCHING,
   };

   [[nodiscard]] LinearSystemScalingMethod linear_system_scaling_method_from_string(const std::string& method);

   // symmetric diagonal scaling s > 0 of a sparse symmetric matrix A given in COO format (any triangle, duplicates are
   // summed): S A S, S = diag(s), is better suited to pivoting. The sparsity pattern is compressed once (initialize),
   // the scaling is computed from the values of the COO entries (compute)
   class SymmetricMatrixScaling {
   public:
      SymmetricMatrixScaling() = default;

      void initialize(size_t dimension, const uno_int* row_indices, const uno_int* column_indices, size_t number_nonzeros,
         int indexing);
      [[nodiscard]] bool is_initialized() const { return this->initialized; }
      [[nodiscard]] bool matches(size_t dimension, size_t number_nonzeros) const {
         return this->dimension == dimension && this->coo_rows.size() == number_nonzeros;
      }

      // compute the scaling vector (of size dimension) from the COO values
      void compute(LinearSystemScalingMethod method, const double* coo_values, std::vector<double>& scaling);

      // scaling factors of the COO entries: s_row * s_column
      void compute_entry_factors(const std::vector<double>& scaling, std::vector<double>& entry_factors) const;

      // parameters of MC77 (ICNTL(7) and CNTL(1)): iterations after the initial pass, and the convergence tolerance on the
      // row norms of the scaled matrix (1 +/- tolerance)
      static constexpr size_t ruiz_max_iterations{10};
      static constexpr double ruiz_tolerance{1e-2};
      static constexpr double min_scaling{1e-100};
      static constexpr double max_scaling{1e100};

   protected:
      bool initialized{false};
      size_t dimension{0};
      // COO entries (0-based) and the index of the compressed entry they contribute to
      std::vector<size_t> coo_rows{};
      std::vector<size_t> coo_columns{};
      std::vector<size_t> coo_to_entry{};
      // compressed (unique) entries of the lower triangle (row >= column)
      std::vector<size_t> entry_rows{};
      std::vector<size_t> entry_columns{};
      std::vector<double> entry_values{};
      // matching-based scaling: assignment (with its workspace), and the column duals of the last matrix (warm start for
      // the next matrix with the same pattern)
      MinimumCostAssignment assignment{};
      std::vector<double> previous_column_duals{};

      void assemble_entry_values(const double* coo_values);
      void compute_ruiz_scaling(std::vector<double>& scaling) const;
      void compute_matching_scaling(std::vector<double>& scaling);
      [[nodiscard]] BipartiteGraph build_matching_graph(std::vector<double>& absolute_values,
         std::vector<double>& log_column_maximum) const;
      void store_column_duals_for_warm_start(const BipartiteGraph& graph);
      static void scale_from_duals(const MinimumCostAssignment& assignment, const std::vector<double>& log_column_maximum,
         const std::vector<size_t>& original_index, std::vector<double>& scaling, std::vector<bool>& scaled_by_matching);
   };
} // namespace

#endif // UNO_SYMMETRICMATRIXSCALING_H
