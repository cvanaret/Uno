// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VIEW_H
#define UNO_VIEW_H

#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>
#include "BLAS.hpp"
#include "symbolic/Inverse.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Subtraction.hpp"
#include "symbolic/symbolic_traits.hpp"
#include "symbolic/Sum.hpp"
#include "symbolic/Transpose.hpp"
#include "symbolic/Triangular.hpp"
#include "symbolic/UnaryNegation.hpp"

namespace uno {
   // constant contiguous array in memory on which BLAS can be called
   template <typename T>
   class View: public SymbolicExpression {
   protected:
      T* pointer;
      size_t view_size;

   public:
      using value_type = std::remove_const_t<T>;

      View(T* pointer, size_t size): pointer(pointer), view_size(size) { }
      View(): pointer(nullptr), view_size(0) { }
      ~View() = default;

      // View<U> -> View<const U> conversion
      template <typename U,
                typename = std::enable_if_t<std::is_convertible_v<U(*)[], T(*)[]>>>
      View(const View<U>& other) noexcept: pointer(other.data()), view_size(other.size()) { }

      View(const View& other) = default;
      View(View&& other) = default;
      View<T>& operator=(const View<T>& other) {
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

      // this returns a pointer to the end of the view
      // careful: this may not be a valid element, do not dereference
      [[nodiscard]] T* end() noexcept {
         return this->pointer + this->size();
      }

      // this returns a pointer to the end of the view
      // careful: this may not be a valid element, do not dereference
      [[nodiscard]] const T* end() const noexcept {
         return this->pointer + this->size();
      }

      [[nodiscard]] T& operator[](size_t index) noexcept {
         assert(index < this->size());
         return this->pointer[index];
      }

      [[nodiscard]] const T& operator[](size_t index) const noexcept {
         assert(index < this->size());
         return this->pointer[index];
      }

      // specialized operation y = x, when the other vector has the member functions data() and size()
      template <typename Vector, typename = std::void_t<decltype(std::declval<const Vector&>().data()),
                                                        decltype(std::declval<const Vector&>().size())>>
      View& operator=(const Vector& other) {
         if (other.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         // dispatch copy function on type of elements
         if constexpr (std::is_same_v<value_type, double>) {
            blas1::copy(this->size(), other.data(), this->data());
         }
         else {
            std::copy(other.data(), other.data() + other.size(), this->data());
         }
         return *this;
      }

      // specialized operation y = -x (note: no BLAS operation available)
      template <typename Vector>
      View& operator=(UnaryNegation<Vector>&& expression) {
         const auto& x = expression.get_expression();
         *this = x;
         this->scale(-1.);
         return *this;
      }

      // specialized operation z = x + a * y (note: no BLAS operation available)
      template <typename Vector1, typename Vector2>
      View& operator=(Sum<Vector1, ScalarMultiple<Vector2>>&& expression) {
         const auto& x = expression.get_left();
         const double a = expression.get_right().get_factor();
         const auto& y = expression.get_right().get_expression();
         *this = y;
         this->scale(a);
         *this += x;
         return *this;
      }

      // specialized operation z = x - a * y (note: no BLAS operation available)
      template <typename Vector1, typename Vector2>
      View& operator=(Subtraction<Vector1, ScalarMultiple<Vector2>>&& expression) {
         const auto& x = expression.get_left();
         const double a = expression.get_right().get_factor();
         const auto& y = expression.get_right().get_expression();
         *this = y;
         this->scale(-a);
         *this += x;
         return *this;
      }

      // specialized operation y = a * x (note: no BLAS operation available)
      template <typename Vector>
      View& operator=(ScalarMultiple<Vector>&& expression) {
         const auto& a = expression.get_factor();
         const auto& x = expression.get_expression();
         if (x.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         *this = x;
         this->scale(a);
         return *this;
      }

      // specialized operation y = x - z
      // note: no BLAS operation available
      template <typename Vector1, typename Vector2>
      View& operator=(Subtraction<Vector1, Vector2>&& expression) {
         const auto& x = expression.get_left();
         const auto& z = expression.get_right();
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
      View& operator=(Multiplication<Matrix, Vector>&& expression) {
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
      View& operator+=(ScalarMultiple<Vector>&& expression) {
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
      View& operator+=(const Vector& other) {
         if (other.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         blas1::add(this->size(), 1., other.data(), this->data());
         return *this;
      }

      // generic operation y -= x
      template <typename Vector>
      View& operator-=(const Vector& other) {
         if (other.size() != this->size()) {
            throw std::invalid_argument("Dimension mismatch between x and y");
         }
         blas1::add(this->size(), -1., other.data(), this->data());
         return *this;
      }

      // x := U⁻ᵀ y with U upper triangular (solve Uᵀ x := y)
      template <typename Matrix, typename Vector>
      View<T>& operator=(Multiplication<Transpose<Inverse<UpperTriangular<Matrix>>>, Vector>&& expression) {
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
      View<T>& operator=(Multiplication<Inverse<UpperTriangular<Matrix>>, Vector>&& expression) {
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
      View& operator=(Multiplication<Transpose<Matrix>, Vector>&& expression) {
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
      View& operator+=(Multiplication<Matrix, Vector>&& expression) {
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
      View& operator-=(Multiplication<Matrix, Vector>&& expression) {
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

      void fill(T value) {
         std::fill(this->data(), this->data() + this->size(), value);
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
      return View{container.data() + start, end - start};
   }

   template <typename T>
   auto view(T* pointer, size_t size) {
      return View{pointer, size};
   }

   template <typename T>
   auto view(T* pointer, size_t start, size_t end) {
      if (start > end) {
         throw std::out_of_range("invalid vector view");
      }
      return View{pointer + start, end - start};
   }

   template <typename T>
   auto view(std::vector<T>& vector) {
      return View{vector.data(), vector.size()};
   }

   template <typename Expression>
   std::ostream& operator<<(std::ostream& stream, const View<Expression>& view) {
      for (size_t index: Range(view.size())) {
         stream << view[index] << " ";
      }
      return stream;
   }

   template <typename Vector1, typename Vector2>
   double dot(const Vector1& x, const Vector2& y) {
      if (x.size() != y.size()) {
         throw std::invalid_argument("Dimension mismatch between x and y");
      }
      return blas1::dot(x.size(), x.data(), y.data());
   }

   inline double dot(const View<double>& x, const double* y) {
      return blas1::dot(x.size(), x.data(), y);
   }
} // namespace

#endif //UNO_VIEW_H