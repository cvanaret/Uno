// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INVERSE_H
#define UNO_INVERSE_H

#include <utility>

namespace uno {
   template <typename Matrix>
   class Inverse {
   public:
      using value_type = typename std::remove_reference_t<Matrix>::value_type;

      explicit Inverse(Matrix&& matrix): matrix(std::forward<Matrix>(matrix)) { }

      [[nodiscard]] const Matrix& get_matrix() const {
         return this->matrix;
      }

   protected:
      const Matrix matrix;
   };

   // free function
   template <typename Matrix>
   inline Inverse<Matrix> inverse(Matrix&& matrix) {
      return Inverse<Matrix>(std::forward<Matrix>(matrix));
   }
} // namespace

#endif // UNO_INVERSE_H