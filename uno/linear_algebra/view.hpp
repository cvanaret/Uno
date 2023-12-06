// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VIEW_H
#define UNO_VIEW_H

// span of an arbitrary container: allocation-free view of a certain length
template <typename ElementType>
class view {
public:
   // compatible with algorithms that query the type of the elements
   using value_type = ElementType;

   view(const std::vector<ElementType>& array, size_t length) noexcept;

   const ElementType& operator[](size_t index) const noexcept;
   [[nodiscard]] size_t size() const noexcept;

   const ElementType* begin() noexcept;
   const ElementType* end() noexcept;

protected:
   const ElementType* array;
   const size_t length;
};

template <typename ElementType>
view<ElementType>::view(const std::vector<ElementType>& array, size_t length) noexcept:
      array(array.data()), length(std::min(length, array.size())) {
}

// preconditions: array != nullptr, i < length
template <typename ElementType>
const ElementType& view<ElementType>::operator[](size_t index) const noexcept {
   return this->array[index];
}

template <typename ElementType>
size_t view<ElementType>::size() const noexcept {
   return this->length;
}

template <typename ElementType>
const ElementType* view<ElementType>::begin() noexcept {
   return this->array;
}

template <typename ElementType>
const ElementType* view<ElementType>::end() noexcept {
   return this->array + this->length;
}

#endif //UNO_VIEW_H