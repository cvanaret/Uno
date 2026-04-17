// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TRIANGULAR_H
#define UNO_TRIANGULAR_H

#include <utility>

namespace uno {
   template <typename Matrix>
   class LowerTriangular {
   public:
      using value_type = typename std::remove_reference_t<Matrix>::value_type;

      explicit LowerTriangular(Matrix&& matrix): matrix(std::forward<Matrix>(matrix)) { }

      [[nodiscard]] const Matrix& get_matrix() const {
         return this->matrix;
      }

   protected:
      const Matrix matrix;
   };

   // free function
   template <typename Matrix>
   inline LowerTriangular<Matrix> lower_triangular(Matrix&& matrix) {
      return LowerTriangular<Matrix>(std::forward<Matrix>(matrix));
   }

   template <typename Matrix>
   class UpperTriangular {
   public:
      using value_type = typename std::remove_reference_t<Matrix>::value_type;

      explicit UpperTriangular(Matrix&& matrix): matrix(std::forward<Matrix>(matrix)) { }

      [[nodiscard]] const Matrix& get_matrix() const {
         return this->matrix;
      }

   protected:
      const Matrix matrix;
   };

   // free function
   template <typename Matrix>
   inline UpperTriangular<Matrix> upper_triangular(Matrix&& matrix) {
      return UpperTriangular<Matrix>(std::forward<Matrix>(matrix));
   }

   template <typename Matrix>
   class LowerUnitTriangular {
   public:
      using value_type = typename std::remove_reference_t<Matrix>::value_type;

      explicit LowerUnitTriangular(Matrix&& matrix): matrix(std::forward<Matrix>(matrix)) { }

      [[nodiscard]] const Matrix& get_matrix() const {
         return this->matrix;
      }

   protected:
      const Matrix matrix;
   };

   // free function
   template <typename Matrix>
   inline LowerUnitTriangular<Matrix> lower_unit_triangular(Matrix&& matrix) {
      return LowerUnitTriangular<Matrix>(std::forward<Matrix>(matrix));
   }

   template <typename Matrix>
   class UpperUnitTriangular {
   public:
      using value_type = typename std::remove_reference_t<Matrix>::value_type;

      explicit UpperUnitTriangular(Matrix&& matrix): matrix(std::forward<Matrix>(matrix)) { }

      [[nodiscard]] const Matrix& get_matrix() const {
         return this->matrix;
      }

   protected:
      const Matrix matrix;
   };

   // free function
   template <typename Matrix>
   inline UpperUnitTriangular<Matrix> upper_unit_triangular(Matrix&& matrix) {
      return UpperUnitTriangular<Matrix>(std::forward<Matrix>(matrix));
   }
} // namespace

#endif // UNO_TRIANGULAR_H