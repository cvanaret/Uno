// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VIEW_H
#define UNO_VIEW_H

// span of an arbitrary container: allocation-free view of a certain length
template <typename T>
class view {
public:
   // compatible with algorithms that query the type of the elements
   using value_type = T;

   view(const std::vector<T>& array, size_t length) noexcept;

   const T& operator[](size_t i) const noexcept;
   [[nodiscard]] size_t size() const noexcept;

   const T* begin() noexcept;
   const T* end() noexcept;

protected:
   const T* array;
   const size_t length;
};

template <typename T>
view<T>::view(const std::vector<T>& array, size_t length) noexcept:
      array(array.data()), length(std::min(length, array.size())) {
}

// preconditions: array != nullptr, i < length
template <typename T>
const T& view<T>::operator[](size_t i) const noexcept {
   return this->array[i];
}

template <typename T>
size_t view<T>::size() const noexcept {
   return this->length;
}

template <typename T>
const T* view<T>::begin() noexcept {
   return this->array;
}

template <typename T>
const T* view<T>::end() noexcept {
   return this->array + this->length;
}

#endif //UNO_VIEW_H