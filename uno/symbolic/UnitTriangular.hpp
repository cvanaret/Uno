// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNITTRIANGULAR_H
#define UNO_UNITTRIANGULAR_H

#include <utility>

namespace uno {
   template <typename Matrix>
   class UnitTriangular {
   public:
      using value_type = typename std::remove_reference_t<Matrix>::value_type;

      explicit UnitTriangular(Matrix&& matrix): matrix(std::forward<Matrix>(matrix)) { }

      [[nodiscard]] const Matrix& get_matrix() const {
         return this->matrix;
      }

   protected:
      const Matrix matrix;
   };

   // free function
   template <typename Matrix>
   inline UnitTriangular<Matrix> unit_triangular(Matrix&& matrix) {
      return UnitTriangular<Matrix>(std::forward<Matrix>(matrix));
   }
} // namespace

#endif // UNO_UNITTRIANGULAR_H