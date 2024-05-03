// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VIEW_H
#define UNO_VIEW_H

// span of an arbitrary container: allocation-free view of a certain length
template <typename Expression>
class view {
public:
   using value_type = typename std::remove_reference_t<Expression>::value_type;

   view(Expression&& array, size_t length) noexcept;

   [[nodiscard]] const typename view::value_type& operator[](size_t index) const noexcept;
   [[nodiscard]] size_t size() const noexcept;

   /*
   const ElementType* begin() noexcept;
   const ElementType* end() noexcept;
    */

protected:
   Expression array;
   const size_t length;
};

template <typename Expression>
view<Expression>::view(Expression&& array, size_t length) noexcept:
      array(std::forward<Expression>(array)), length(std::min(length, array.size())) {
}

// preconditions: array != nullptr, i < length
template <typename Expression>
const typename view<Expression>::value_type& view<Expression>::operator[](size_t index) const noexcept {
   return this->array[index];
}

template <typename Expression>
size_t view<Expression>::size() const noexcept {
   return this->length;
}

/*
template <typename ElementType>
const ElementType* view<ElementType>::begin() noexcept {
   return this->array;
}

template <typename ElementType>
const ElementType* view<ElementType>::end() noexcept {
   return this->array + this->length;
}
*/

#endif //UNO_VIEW_H