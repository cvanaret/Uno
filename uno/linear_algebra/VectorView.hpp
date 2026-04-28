// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTORVIEW_H
#define UNO_VECTORVIEW_H

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include "linear_algebra/BLASVector.hpp"

namespace uno {
   // span of an arbitrary container: allocation-free view of a certain length
   template <typename Vector>
   class VectorView: public BLASVector<typename std::remove_reference_t<Vector>::value_type> {
   public:
      using value_type = typename std::remove_reference_t<Vector>::value_type;

      VectorView(const Vector& vector, size_t start, size_t end):
            BLASVector<value_type>(), vector(vector), start(start), end(std::min(end, this->vector.size())) {
         if (this->end < this->start) {
            throw std::runtime_error("The view ends before its starting point");
         }
      }
      VectorView(const Vector&&) = delete;

      [[nodiscard]] size_t size() const noexcept override {
         return this->end - this->start;
      }

      [[nodiscard]] const value_type* data() const noexcept override {
         return this->vector.data() + this->start;
      }

      [[nodiscard]] const value_type& operator[](size_t index) const noexcept override {
         return this->vector[this->start + index];
      }

   protected:
      const Vector& vector;
      const size_t start;
      const size_t end;
   };

   // span of an arbitrary container: allocation-free view of a certain length
   template <typename Vector>
   class MutableVectorView: public MutableBLASVector<typename std::remove_reference_t<Vector>::value_type> {
   public:
      using value_type = typename std::remove_reference_t<Vector>::value_type;

      MutableVectorView(Vector& vector, size_t start, size_t end):
            MutableBLASVector<value_type>(), vector(vector), start(start), end(std::min(end, this->vector.size())) {
         if (this->end < this->start) {
            throw std::runtime_error("The view ends before its starting point");
         }
      }
      MutableVectorView(const Vector&&) = delete;

      // operators
      using MutableBLASVector<value_type>::operator=;
      using MutableBLASVector<value_type>::operator+=;
      using MutableBLASVector<value_type>::operator-=;

      [[nodiscard]] size_t size() const noexcept override {
         return this->end - this->start;
      }

      [[nodiscard]] value_type* data() noexcept override {
         return this->vector.data() + this->start;
      }

      [[nodiscard]] const value_type* data() const noexcept override {
         return this->vector.data() + this->start;
      }

      [[nodiscard]] value_type& operator[](size_t index) noexcept override {
         return this->vector[this->start + index];
      }

      [[nodiscard]] const value_type& operator[](size_t index) const noexcept override {
         return this->vector[this->start + index];
      }

   protected:
      Vector& vector;
      const size_t start;
      const size_t end;
   };

   // span of a C array: allocation-free view of a certain length
   template <typename T>
   class MutableArrayView: public MutableBLASVector<T> {
   public:
      MutableArrayView(T* array, size_t start, size_t end):
            array(array), start(start), end(end) {
         if (this->end < this->start) {
            throw std::runtime_error("The view ends before its starting point");
         }
      }

      // operators
      using MutableBLASVector<T>::operator=;
      using MutableBLASVector<T>::operator+=;
      using MutableBLASVector<T>::operator-=;

      [[nodiscard]] size_t size() const noexcept override {
         return this->end - this->start;
      }

      [[nodiscard]] T* data() noexcept override {
         return this->array + this->start;
      }

      [[nodiscard]] const T* data() const noexcept override {
         return this->array + this->start;
      }

      [[nodiscard]] T& operator[](size_t index) noexcept override {
         return this->array[this->start + index];
      }

      [[nodiscard]] const T& operator[](size_t index) const noexcept override {
         return this->array[this->start + index];
      }

   protected:
      T* array;
      const size_t start;
      const size_t end;
   };

   // free functions
   template <typename Expression>
   auto view(const Expression& expression, size_t start, size_t end) {
      return VectorView{expression, start, end};
   }

   template <typename Expression>
   auto view(Expression& expression, size_t start, size_t end) {
      return MutableVectorView{expression, start, end};
   }

   template <typename T>
   auto view(T* array, size_t start, size_t end) {
      return MutableArrayView{array, start, end};
   }

   template <typename Expression>
   std::ostream& operator<<(std::ostream& stream, const VectorView<Expression>& view) {
      for (size_t index: Range(view.size())) {
         stream << view[index] << " ";
      }
      return stream;
   }
} // namespace

#endif //UNO_VECTORVIEW_H