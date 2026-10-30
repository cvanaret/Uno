// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICMATRIXSCALING_H
#define UNO_SYMMETRICMATRIXSCALING_H

#include <cstddef>
#include <string>
#include <vector>
#include "../interfaces/C/uno_int.h"

namespace uno {
   enum class LinearSystemScalingMethod {
      NONE,
      RUIZ, // symmetric Ruiz infinity-norm equilibration (Ruiz 2001; Knight, Ruiz & Uçar 2014)
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

      // parameters
      static constexpr size_t ruiz_max_iterations{10};
      static constexpr double ruiz_tolerance{1e-2}; // stop when all the row norms are in [1 - tol, 1 + tol]
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

      void assemble_entry_values(const double* coo_values);
      void compute_ruiz_scaling(std::vector<double>& scaling) const;
   };
} // namespace

#endif // UNO_SYMMETRICMATRIXSCALING_H
