// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BLASVECTOR_H
#define UNO_BLASVECTOR_H

#include <algorithm>
#include <cstddef>
#include <cassert>
#include "BLAS.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Subtraction.hpp"
#include "symbolic/Transpose.hpp"

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
   };

   // non-constant contiguous array in memory on which BLAS can be called
   template <typename T>
   class MutableBLASVector: public BLASVector<T> {
   public:
      MutableBLASVector() = default;

      // specialized operation y = x, when the other vector has the member function data()
      template <typename Vector, decltype(Vector{}.data()) = true>
      auto operator=(const Vector& other) {
         assert(other.size() == this->size());
         blas1::copy(this->size(), other.data(), this->data());
         return *this;
      }

      // generic operation y = expression
      template <typename Expression>
      MutableBLASVector& operator=(const Expression& expression) {
         static_assert(std::is_same_v<typename Expression::value_type, T>);
         assert(expression.size() == this->size());
         for (size_t index: Range(expression.size())) {
            this->operator[](index) = expression[index];
         }
         return *this;
      }

      // specialized operation y = a * x
      // note: no BLAS operation available
      template <typename Vector>
      MutableBLASVector& operator=(ScalarMultiple<Vector>&& expression) {
         const auto& a = expression.get_factor();
         const auto& x = expression.get_expression();
         assert(x.size() == this->size());
         for (size_t index: Range(this->size())) {
            this->operator[](index) = a * x[index];
         }
         return *this;
      }

      // specialized operation y = x - z
      // note: no BLAS operation available
      template <typename Vector>
      MutableBLASVector& operator=(Subtraction<Vector, Vector>&& expression) {
         *this = expression.get_expression1();
         *this -= expression.get_expression2();
         return *this;
      }

      // specialized operation y += a * x
      template <typename Vector>
      MutableBLASVector& operator+=(ScalarMultiple<Vector>&& expression) {
         const auto& a = expression.get_factor();
         const auto& x = expression.get_expression();
         assert(x.size() == this->size());
         blas1::add(this->size(), a, x.data(), this->data());
         return *this;
      }

      // generic operation y += x
      template <typename Vector>
      MutableBLASVector& operator+=(const Vector& other) {
         assert(other.size() == this->size());
         blas1::add(this->size(), 1., other.data(), this->data());
         return *this;
      }

      template <typename Vector>
      MutableBLASVector& operator-=(const Vector& other) {
         assert(other.size() == this->size());
         blas1::add(this->size(), -1., other.data(), this->data());
         return *this;
      }

      // y := A^T x
      template <typename Matrix, typename Vector>
      MutableBLASVector& operator=(Multiplication<Transpose<Matrix>, Vector>&& expression) {
         const auto& A = expression.get_left().get_matrix();
         const auto& x = expression.get_right();
         assert(A.number_rows == x.size());
         blas2::matrix_vector_product('T', A.number_rows, A.number_columns, 1., A.data(), A.leading_dimension, x.data(),
            0., this->data());
         return *this;
      }

      // y -= Ax (or y = -A x + y)
      template <typename Matrix, typename Vector>
      MutableBLASVector& operator-=(Multiplication<Matrix, Vector>&& expression) {
         const auto& A = expression.get_left();
         const auto& x = expression.get_right();
         assert(A.number_columns == x.size());
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

   inline double dot(const BLASVector<double>& x, const BLASVector<double>& y) {
      return blas1::dot(std::min(x.size(), y.size()), x.data(), y.data());
   }

   inline double dot(const BLASVector<double>& x, const double* y) {
      return blas1::dot(x.size(), x.data(), y);
   }
} // namespace

#endif //UNO_BLASVECTOR_H