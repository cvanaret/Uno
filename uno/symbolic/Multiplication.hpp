// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MULTIPLICATION_H
#define UNO_MULTIPLICATION_H

#include <utility>

namespace uno {
   // stores the expression (matrix1 * matrix2) symbolically
   template <typename Matrix1, typename Matrix2>
   class Multiplication {
   public:
      using value_type = typename std::remove_reference_t<Matrix1>::value_type;

      Multiplication(Matrix1&& expression1, Matrix2&& matrix2):
            matrix1(std::forward<Matrix1>(expression1)), matrix2(std::forward<Matrix2>(matrix2)) {
      }

      [[nodiscard]] const Matrix1& get_left() const {
         return this->matrix1;
      }

      [[nodiscard]] const Matrix2& get_right() const {
         return this->matrix2;
      }

   protected:
      const Matrix1 matrix1;
      const Matrix2 matrix2;
   };

   // free function
   template <typename Matrix1, typename Matrix2,
         std::enable_if_t<std::is_same_v<typename std::remove_reference_t<Matrix1>::value_type,
                                         typename std::remove_reference_t<Matrix2>::value_type>, int> = 0,
      // Matrix1 and Matrix2 are both not arithmetic types
      std::enable_if_t<!std::is_arithmetic_v<Matrix1>, int> = 0,
      std::enable_if_t<!std::is_arithmetic_v<Matrix2>, int> = 0>
   inline Multiplication<Matrix1, Matrix2> operator*(Matrix1&& matrix1, Matrix2&& matrix2) {
      return {std::forward<Matrix1>(matrix1), std::forward<Matrix2>(matrix2)};
   }
} // namespace

#endif // UNO_MULTIPLICATION_H