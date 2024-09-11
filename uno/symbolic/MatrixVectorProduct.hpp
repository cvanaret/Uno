// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MATRIXVECTORPRODUCT_H
#define UNO_MATRIXVECTORPRODUCT_H

#include "linear_algebra/SparseVector.hpp"

namespace uno {
   // symbolic matrix-vector product
   template <typename Matrix, typename Vector>
   class MatrixVectorProduct {
   public:
      using value_type = typename std::remove_reference_t<Vector>::value_type;

      MatrixVectorProduct(Matrix&& matrix, Vector&& vector): matrix(std::forward<Matrix>(matrix)), vector(std::forward<Vector>(vector)) { }

      [[nodiscard]] constexpr size_t size() const { return this->vector.size(); }

      // product computed using row-major matrix
      [[nodiscard]] typename MatrixVectorProduct::value_type operator[](size_t row_index) const {
         return dot(this->vector, this->matrix[row_index]);
      }

   protected:
      Matrix matrix;
      Vector vector;
   };

   // free function
   template <typename Matrix, typename Vector,
         typename std::enable_if<std::is_same_v<typename std::remove_reference_t<Matrix>::value_type,
                                                typename std::remove_reference_t<Vector>::value_type>, int>::type = 0>
   inline MatrixVectorProduct<Matrix, Vector> operator*(Matrix&& matrix, Vector&& vector) {
      return {std::forward<Matrix>(matrix), std::forward<Vector>(vector)};
   }
} // namespace

#endif // UNO_MATRIXVECTORPRODUCT_H
