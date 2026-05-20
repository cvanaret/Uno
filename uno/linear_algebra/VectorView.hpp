// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTORVIEW_H
#define UNO_VECTORVIEW_H

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
   class VectorView {
   protected:
      T* pointer;
      const size_t view_size;

   public:
      using value_type = std::remove_const_t<T>;

      VectorView(T* pointer, size_t size): pointer(pointer), view_size(size) { }
      ~VectorView() = default;
      VectorView(const VectorView& other) = default;
      VectorView(VectorView&& other) = default;
      VectorView<T>& operator=(const VectorView<T>& other) {
         if (&other != this) {
            blas1::copy(this->size(), other.data(), this->data());
         }
         return *this;
      }

      [[nodiscard]] size_t size() const noexcept {
         return this->view_size;
      }

      [[nodiscard]] bool empty() const {
         return (this->size() == 0);
      }

      [[nodiscard]] T* data() noexcept {
         return this->pointer;
      }

      [[nodiscard]] const T* data() const noexcept {
         return this->pointer;
      }

      [[nodiscard]] T& operator[](size_t index) noexcept {
         return this->pointer[index];
      }

      [[nodiscard]] const T& operator[](size_t index) const noexcept {
         return this->pointer[index];
      }

      // specialized operation y = x, when the other vector has the member function data()
      template <typename Vector, decltype(Vector{}.data()) = true>
      VectorView& operator=(const Vector& other) {
         if (other.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         blas1::copy(this->size(), other.data(), this->data());
         return *this;
      }

      // generic operation y = expression
      template <typename Expression>
      VectorView& operator=(const Expression& expression) {
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
      VectorView& operator=(ScalarMultiple<Vector>&& expression) {
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
      VectorView& operator=(Subtraction<Vector, Vector>&& expression) {
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
      VectorView& operator=(Multiplication<Matrix, Vector>&& expression) {
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
      VectorView& operator+=(ScalarMultiple<Vector>&& expression) {
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
      VectorView& operator+=(const Vector& other) {
         if (other.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         blas1::add(this->size(), 1., other.data(), this->data());
         return *this;
      }

      // generic operation y -= x
      template <typename Vector>
      VectorView& operator-=(const Vector& other) {
         if (other.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         blas1::add(this->size(), -1., other.data(), this->data());
         return *this;
      }

      // x := U⁻ᵀ y with U upper triangular (solve Uᵀ x := y)
      template <typename Matrix, typename Vector>
      VectorView<T>& operator=(Multiplication<Transpose<Inverse<UpperTriangular<Matrix>>>, Vector>&& expression) {
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
      VectorView<T>& operator=(Multiplication<Inverse<UpperTriangular<Matrix>>, Vector>&& expression) {
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
      VectorView& operator=(Multiplication<Transpose<Matrix>, Vector>&& expression) {
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
      VectorView& operator+=(Multiplication<Matrix, Vector>&& expression) {
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
      VectorView& operator-=(Multiplication<Matrix, Vector>&& expression) {
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

      void scale(T factor) {
         blas1::scale(this->size(), factor, this->data());
      }

      void operator*=(T factor) {
         this->scale(factor);
      }
   };

   // free functions
   template <typename Container>
   auto view(Container& container, size_t start, size_t end) {
      if (start > end || end > container.size()) {
         throw std::out_of_range("invalid vector view");
      }
      return VectorView{container.data() + start, end - start};
   }

   template <typename T>
   auto view(T* pointer, size_t size) {
      return VectorView{pointer, size};
   }

   template <typename T>
   auto view(T* pointer, size_t start, size_t end) {
      if (start > end) {
         throw std::out_of_range("invalid vector view");
      }
      return VectorView{pointer + start, end - start};
   }

   template <typename Expression>
   std::ostream& operator<<(std::ostream& stream, const VectorView<Expression>& view) {
      for (size_t index: Range(view.size())) {
         stream << view[index] << " ";
      }
      return stream;
   }

   template <typename V1, typename V2>
   double dot(const V1& x, const V2& y) {
      if (x.size() != y.size()) {
         throw std::invalid_argument("Dimension mismatch between x and y");
      }
      return blas1::dot(x.size(), x.data(), y.data());
   }

   inline double dot(const VectorView<double>& x, const double* y) {
      return blas1::dot(x.size(), x.data(), y);
   }
} // namespace

#endif //UNO_VECTORVIEW_H