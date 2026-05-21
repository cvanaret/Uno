// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TRANSPOSE_H
#define UNO_TRANSPOSE_H

#include "symbolic_traits.hpp"

namespace uno {
   template <typename Matrix>
   class Transpose {
   public:
      using value_type = typename std::remove_reference_t<Matrix>::value_type;

      explicit Transpose(Matrix&& matrix): matrix(std::forward<Matrix>(matrix)) { }

      [[nodiscard]] constexpr decltype(auto) get_matrix() const noexcept {
         return this->matrix;
      }

   protected:
      storage_t<Matrix> matrix;
   };

   // free function
   template <typename Matrix>
   Transpose<Matrix> transpose(Matrix&& matrix) {
      return Transpose<Matrix>(std::forward<Matrix>(matrix));
   }
} // namespace

#endif // UNO_TRANSPOSE_H