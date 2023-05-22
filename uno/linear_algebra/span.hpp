// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SPAN_H
#define UNO_SPAN_H

// span of an arbitrary container: allocation-free view of a certain length
template <typename T>
class span {
public:
   // compatible with algorithms that query the type of the elements
   using value_type = T;

   span(const std::vector<T>& array, size_t length) noexcept;

   const T& operator[](size_t i) const noexcept;
   [[nodiscard]] size_t size() const noexcept;

   const T* begin() noexcept;
   const T* end() noexcept;

protected:
   const T* array;
   const size_t length;
};

template <typename T>
span<T>::span(const std::vector<T>& array, size_t length) noexcept:
      array(array.data()), length(std::min(length, array.size())) {
}

// preconditions: array != nullptr, i < length
template <typename T>
const T& span<T>::operator[](size_t i) const noexcept {
   return this->array[i];
}

template <typename T>
size_t span<T>::size() const noexcept {
   return this->length;
}

template <typename T>
const T* span<T>::begin() noexcept {
   return this->array;
}

template <typename T>
const T* span<T>::end() noexcept {
   return this->array + this->length;
}

#endif //UNO_SPAN_H