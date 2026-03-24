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
         assert(other.size() <= this->size() && "The other vector is larger than the current vector");
         const int size = static_cast<int>(this->size());
         constexpr int increment = 1;
         BLAS_copy_vector(&size, other.data(), &increment, this->data(), &increment);
         return *this;
      }

      // generic operation y = x
      template <typename Vector>
      MutableBLASVector& operator=(const Vector& other) {
         static_assert(std::is_same_v<typename Vector::value_type, T>);
         assert(other.size() <= this->size() && "The expression is larger than the current vector");
         for (size_t index: Range(other.size())) {
            this->operator[](index) = other[index];
         }
         return *this;
      }

      // specialized operation y = a * x
      // note: no BLAS operation available
      template <typename Vector>
      MutableBLASVector& operator=(ScalarMultiple<Vector>&& other) {
         assert(other.get_expression().size() <= this->size() && "The other vector is larger than the current vector");
         for (size_t index: Range(this->size())) {
            this->operator[](index) = other.get_factor() * other.get_expression()[index];
         }
         return *this;
      }

      // specialized operation y = x - z
      // note: no BLAS operation available
      template <typename Vector>
      MutableBLASVector& operator=(Subtraction<Vector, Vector>&& other) {
         *this = other.get_expression1();
         *this -= other.get_expression2();
         return *this;
      }

      // specialized operation y += a * x
      template <typename Vector>
      MutableBLASVector& operator+=(ScalarMultiple<Vector>&& other) {
         assert(other.size() <= this->size() && "The other vector is larger than the current vector");
         const int size = static_cast<int>(this->size());
         const double factor = other.get_factor();
         constexpr int increment = 1;
         BLAS_add_vectors(&size, &factor, other.get_expression().data(), &increment, this->data(), &increment);
         return *this;
      }

      // generic operation y += x
      template <typename Vector>
      MutableBLASVector& operator+=(const Vector& other) {
         assert(other.size() <= this->size() && "The other vector is larger than the current vector");
         const int size = static_cast<int>(this->size());
         constexpr double factor = 1.;
         constexpr int increment = 1;
         BLAS_add_vectors(&size, &factor, other.data(), &increment, this->data(), &increment);
         return *this;
      }

      template <typename Vector>
      MutableBLASVector& operator-=(const Vector& other) {
         assert(other.size() <= this->size() && "The other vector is larger than the current vector");
         const int size = static_cast<int>(this->size());
         constexpr double factor = -1.;
         constexpr int increment = 1;
         BLAS_add_vectors(&size, &factor, other.data(), &increment, this->data(), &increment);
         return *this;
      }

      // y := A^T x
      template <typename Matrix, typename Vector>
      MutableBLASVector& operator=(Multiplication<Transpose<Matrix>, Vector>&& expression) {
         const auto& A = expression.get_left().get_matrix();
         const auto& x = expression.get_right();
         assert(A.number_rows == x.size());
         constexpr char trans = 'T';
         const int m = A.number_rows;
         const int n = A.number_columns;
         constexpr double alpha = 1.;
         const int lda = A.leading_dimension;
         constexpr int increment = 1;
         constexpr double beta = 0.;
         BLAS_matrix_vector_product(&trans, &m, &n, &alpha, A.data(), &lda, x.data(), &increment, &beta, this->data(), &increment);
         return *this;
      }

      // y -= Ax (or y = -A x + y)
      template <typename Matrix, typename Vector>
      MutableBLASVector& operator-=(Multiplication<Matrix, Vector>&& expression) {
         const auto& A = expression.get_left();
         const auto& x = expression.get_right();
         assert(A.number_columns == x.size());
         constexpr char trans = 'N';
         const int m = A.number_rows;
         const int n = A.number_columns;
         constexpr double alpha = -1.;
         const int lda = A.leading_dimension;
         constexpr int increment = 1;
         constexpr double beta = 1.;
         BLAS_matrix_vector_product(&trans, &m, &n, &alpha, A.data(), &lda, x.data(), &increment, &beta, this->data(), &increment);
         return *this;
      }

      [[nodiscard]] virtual T* data() = 0;
      [[nodiscard]] virtual T& operator[](size_t index) = 0;

      void scale(T factor) {
         const int size = static_cast<int>(this->size());
         constexpr int increment = 1;
         BLAS_scale_vector(&size, &factor, this->data(), &increment);
      }

      void operator*=(T factor) {
         this->scale(factor);
      }
   };

   inline double dot(const BLASVector<double>& x, const BLASVector<double>& y) {
      const int size = static_cast<int>(std::min(x.size(), y.size()));
      constexpr int increment = 1;
      return BLAS_dot_product(&size, x.data(), &increment, y.data(), &increment);
   }

   inline double dot(const BLASVector<double>& x, const double* y) {
      const int size = static_cast<int>(x.size());
      constexpr int increment = 1;
      return BLAS_dot_product(&size, x.data(), &increment, y, &increment);
   }
} // namespace

#endif //UNO_BLASVECTOR_H