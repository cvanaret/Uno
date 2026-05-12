// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTORVIEW_H
#define UNO_VECTORVIEW_H

#include <cstddef>
#include <stdexcept>
#include "linear_algebra/BLASVector.hpp"

namespace uno {
   // span of a C pointer: allocation-free view of a certain length
   template <typename T>
   class VectorView: public BLASVector<T> {
   public:
      using value_type = T;

      VectorView(T* pointer, size_t size) noexcept: pointer(pointer), size_(size) {
      }

      // operators
      using BLASVector<T>::operator=;
      using BLASVector<T>::operator+=;
      using BLASVector<T>::operator-=;

      [[nodiscard]] size_t size() const noexcept override {
         return this->size_;
      }

      [[nodiscard]] T* data() noexcept override {
         return this->pointer;
      }

      [[nodiscard]] const T* data() const noexcept override {
         return this->pointer;
      }

      [[nodiscard]] T& operator[](size_t index) noexcept override {
         return this->pointer[index];
      }

      [[nodiscard]] const T& operator[](size_t index) const noexcept override {
         return this->pointer[index];
      }

   protected:
      T* pointer;
      const size_t size_;
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
   auto view(T* pointer, size_t start, size_t end) {
      return VectorView{pointer + start, end - start};
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