// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BLASVECTOR_H
#define UNO_BLASVECTOR_H

#include <cstddef>
#include "BLAS.hpp"
#include "symbolic/Inverse.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Subtraction.hpp"
#include "symbolic/Transpose.hpp"
#include "symbolic/Triangular.hpp"

namespace uno {
   // constant contiguous array in memory on which BLAS can be called
   template <typename T>
   class BLASVector {
   public:
      using value_type = T;

      BLASVector() = default;
      virtual ~BLASVector() = default;

      [[nodiscard]] virtual size_t size() const = 0;
      [[nodiscard]] bool empty() const {
         return (this->size() == 0);
      }
      [[nodiscard]] virtual const T* data() const = 0;
      [[nodiscard]] virtual const T& operator[](size_t index) const = 0;

      // specialized operation y = x, when the other vector has the member function data()
      template <typename Vector, decltype(Vector{}.data()) = true>
      auto operator=(const Vector& other) {
         if (other.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         blas1::copy(this->size(), other.data(), this->data());
         return *this;
      }

      // generic operation y = expression
      template <typename Expression>
      BLASVector& operator=(const Expression& expression) {
         // static_assert(std::is_same_v<typename Expression::value_type, T>);
         if (expression.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between expression and y");
         }
         for (size_t index: Range(expression.size())) {
            this->operator[](index) = expression[index];
         }
         return *this;
      }

      // specialized operation y = a * x (note: no BLAS operation available)
      template <typename Vector>
      BLASVector& operator=(ScalarMultiple<Vector>&& expression) {
         const auto& a = expression.get_factor();
         const auto& x = expression.get_expression();
         if (x.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         for (size_t index: Range(this->size())) {
            this->operator[](index) = a * x[index];
         }
         return *this;
      }

      // specialized operation y = x - z
      // note: no BLAS operation available
      template <typename Vector>
      BLASVector& operator=(Subtraction<Vector, Vector>&& expression) {
         const auto& x = expression.get_expression1();
         const auto& z = expression.get_expression2();
         if (x.size() != z.size()) {
            throw std::invalid_argument("Dimension mismatch between x and z");
         }
         if (x.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         *this = x;
         *this -= z;
         return *this;
      }

      // y = Ax (or y = A x + 0 * y)
      template <typename Matrix, typename Vector>
      BLASVector& operator=(Multiplication<Matrix, Vector>&& expression) {
         const auto& x = expression.get_right();
         if (x.data() == this->data()) {
            throw std::invalid_argument("BLASVector::operator= cannot be called with x == y");
         }
         const auto& A = expression.get_left();
         if (A.number_columns != x.size()) {
            throw std::invalid_argument("Dimension mismatch between A and x");
         }
         if (A.number_rows != this->size()) {
            throw std::invalid_argument("Dimension mismatch between A and y");
         }
         blas2::matrix_vector_product('N', A.number_rows, A.number_columns, 1., A.data(), A.leading_dimension, x.data(),
            0., this->data());
         return *this;
      }

      // specialized operation y += a * x
      template <typename Vector>
      BLASVector& operator+=(ScalarMultiple<Vector>&& expression) {
         const auto& a = expression.get_factor();
         const auto& x = expression.get_expression();
         if (x.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         blas1::add(this->size(), a, x.data(), this->data());
         return *this;
      }

      // generic operation y += x
      template <typename Vector>
      BLASVector& operator+=(const Vector& other) {
         if (other.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         blas1::add(this->size(), 1., other.data(), this->data());
         return *this;
      }

      // generic operation y -= x
      template <typename Vector>
      BLASVector& operator-=(const Vector& other) {
         if (other.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         blas1::add(this->size(), -1., other.data(), this->data());
         return *this;
      }

      // x := U⁻ᵀ y with U upper triangular (solve Uᵀ x := y)
      template <typename Matrix, typename Vector>
      BLASVector<T>& operator=(Multiplication<Transpose<Inverse<UpperTriangular<Matrix>>>, Vector>&& expression) {
         const auto& U = expression.get_left().get_matrix().get_matrix().get_matrix();
         const auto& y = expression.get_right();
         if (this->size() != U.number_columns) {
            throw std::runtime_error("Dimension mismatch in BLASVector::operator=");
         }
         // copy the RHS into this
         blas1::copy(this->size(), y.data(), this->data());
         blas3::triangular_back_solve('L', 'U', 'T' /* transpose */, 'N' /* non-unit */, y.size(), 1, 1., U.data(),
            U.leading_dimension, this->data(), this->size());
         return *this;
      }

      // x := U⁻¹ y with U upper triangular (solve U x := y)
      template <typename Matrix, typename Vector>
      BLASVector<T>& operator=(Multiplication<Inverse<UpperTriangular<Matrix>>, Vector>&& expression) {
         const auto& U = expression.get_left().get_matrix().get_matrix();
         const auto& y = expression.get_right();
         if (this->size() != U.number_columns) {
            throw std::runtime_error("Dimension mismatch in BLASVector::operator=");
         }
         // copy the RHS into this
         blas1::copy(this->size(), y.data(), this->data());
         blas3::triangular_back_solve('L', 'U', 'N' /* no transpose */, 'N' /* non-unit */, y.size(), 1, 1., U.data(),
            U.leading_dimension, this->data(), this->size());
         return *this;
      }

      // y := A^T x
      template <typename Matrix, typename Vector>
      BLASVector& operator=(Multiplication<Transpose<Matrix>, Vector>&& expression) {
         const auto& A = expression.get_left().get_matrix();
         const auto& x = expression.get_right();
         if (A.number_rows != x.size()) {
            throw std::invalid_argument("Dimension mismatch between A and x");
         }
         if (A.number_columns != this->size()) {
            throw std::invalid_argument("Dimension mismatch between A and y");
         }
         blas2::matrix_vector_product('T', A.number_rows, A.number_columns, 1., A.data(), A.leading_dimension, x.data(),
            0., this->data());
         return *this;
      }

      // y += Ax (or y = A x + y)
      template <typename Matrix, typename Vector>
      BLASVector& operator+=(Multiplication<Matrix, Vector>&& expression) {
         const auto& A = expression.get_left();
         const auto& x = expression.get_right();
         if (A.number_columns != x.size()) {
            throw std::invalid_argument("Dimension mismatch between A and x");
         }
         if (A.number_rows != this->size()) {
            throw std::invalid_argument("Dimension mismatch between A and y");
         }
         blas2::matrix_vector_product('N', A.number_rows, A.number_columns, 1., A.data(), A.leading_dimension, x.data(),
            1., this->data());
         return *this;
      }

      // y -= Ax (or y = -A x + y)
      template <typename Matrix, typename Vector>
      BLASVector& operator-=(Multiplication<Matrix, Vector>&& expression) {
         const auto& A = expression.get_left();
         const auto& x = expression.get_right();
         if (A.number_columns != x.size()) {
            throw std::invalid_argument("Dimension mismatch between A and x");
         }
         if (A.number_rows != this->size()) {
            throw std::invalid_argument("Dimension mismatch between A and y");
         }
         blas2::matrix_vector_product('N', A.number_rows, A.number_columns, -1., A.data(), A.leading_dimension, x.data(),
            1., this->data());
         return *this;
      }

      [[nodiscard]] virtual T* data() = 0;

      [[nodiscard]] virtual T& operator[](size_t index) = 0;

      void scale(T factor) {
         blas1::scale(this->size(), factor, this->data());
      }

      void operator*=(T factor) {
         this->scale(factor);
      }
   };

   template <typename T, typename U>
   double dot(const BLASVector<T>& x, const BLASVector<U>& y) {
      if (x.size() != y.size()) {
         throw std::invalid_argument("Dimension mismatch between x and y");
      }
      return blas1::dot(x.size(), x.data(), y.data());
   }

   inline double dot(const BLASVector<double>& x, const double* y) {
      return blas1::dot(x.size(), x.data(), y);
   }
} // namespace

#endif //UNO_BLASVECTOR_H