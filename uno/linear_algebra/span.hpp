// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SPAN_H
#define UNO_SPAN_H

// span of an arbitrary container: allocation-free view of a certain length
template <typename ARRAY>
class span {
public:
   // compatible with algorithms that query the type of the elements
   using value_type = typename ARRAY::value_type;

   span(const ARRAY& array, size_t length) noexcept;
   const value_type& operator[](size_t i) const noexcept;
   [[nodiscard]] size_t size() const noexcept;
   const ARRAY& begin() noexcept;
   const ARRAY& end() noexcept;

protected:
   const ARRAY& array;
   size_t length;
};

template <typename ARRAY>
span<ARRAY>::span(const ARRAY& array, size_t length) noexcept:
      array(array), length(std::min(length, array.size())) {
}

template <typename ARRAY>
const typename span<ARRAY>::value_type& span<ARRAY>::operator[](size_t i) const noexcept {
   return this->array[i];
}

template <typename ARRAY>
size_t span<ARRAY>::size() const noexcept {
   return this->length;
}

template <typename ARRAY>
const ARRAY& span<ARRAY>::begin() noexcept {
   return *this;
}

template <typename ARRAY>
const ARRAY& span<ARRAY>::end() noexcept {
   return *this;
}

#endif //UNO_SPAN_H