// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAPACK_EXTENSION_H
#define UNO_LAPACK_EXTENSION_H

#include <cmath>
#include "LAPACK.hpp"
#include "VectorView.hpp"

namespace uno {
   // LDL' (aka signed Cholesky) factorization without pivoting for symmetric indefinite matrices
   // Returns true upon success (all pivots are larger in norm than zero_pivot_tolerance), false upon failure
   // Partition A (n×n, column-major with leading dimension lda, symmetric) as
   //
   //   A = [ a11  a21ᵀ ]
   //       [ a21  A22  ]
   //
   // The LDL' decomposition factorises the leading 1×1 block first:
   //
   //   d11  = a11
   //   l21  = a21 / d11
   //   A22' = A22 - d11 · l21 · l21ᵀ          (Schur complement)
   //
   // then recurses on the (n-1)×(n-1) trailing Schur complement A22'.
   // Only the lower triangle of A is read and updated throughout.
   inline bool ldlt_nopiv_lvl2_rightlooking(double* A, size_t n, size_t lda, double zero_pivot_tolerance) noexcept {
      // base case
      if (n == 0) {
         return true;
      }

      // step 1: extract pivot d11 = A(0,0)
      const double d11 = A[0];
      // near zero pivot: failure
      if (std::abs(d11) <= zero_pivot_tolerance) {
         return false;
      }

      // step 2: compute column of L:  l21 <- a21 / d11
      auto l21 = view(A, 1, n - 1);
      l21.scale(1./d11);

      // step 3: rank-1 update of the trailing sub-matrix (lower triangle):
      //   A22 <- A22 - d11 · l21 · l21ᵀ
      auto A22 = view(A, 1 + lda, lda*n - 1);
      blas3::symmetric_rank_1_update('L', l21.size(), -d11, l21.data(), A22.data(), lda);

      // step 4: recurse on the (n-1)×(n-1) trailing sub-matrix
      return ldlt_nopiv_lvl2_rightlooking(A22.data(), n - 1, lda, zero_pivot_tolerance);
   }
} // namespace

#endif // UNO_LAPACK_EXTENSION_H