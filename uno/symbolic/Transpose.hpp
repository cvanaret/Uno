// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TRANSPOSE_H
#define UNO_TRANSPOSE_H

#include <cstddef>
#include <utility>

namespace uno {
   template <typename Matrix>
   class Transpose {
   public:
      using value_type = typename std::remove_reference_t<Matrix>::value_type;

      explicit Transpose(Matrix&& matrix): matrix(std::forward<Matrix>(matrix)) { }

      [[nodiscard]] constexpr size_t size() const {
         return this->matrix.size();
      }

      [[nodiscard]] const Matrix& get_matrix() const {
         return this->matrix;
      }

   protected:
      const Matrix matrix;
   };

   // free function
   template <typename Matrix>
   inline Transpose<Matrix> transpose(Matrix&& matrix) {
      return Transpose<Matrix>(std::forward<Matrix>(matrix));
   }
} // namespace

#endif // UNO_TRANSPOSE_H
