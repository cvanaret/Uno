// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTORVIEW_H
#define UNO_VECTORVIEW_H

#include <vector>

template <typename VectorType>
class VectorView {
public:
   using value_type = typename VectorType::value_type;

   // end index is excluded
   VectorView(VectorType& vector, size_t start_index, size_t end_index);
   value_type& operator[](size_t index);
   const value_type& operator[](size_t index) const;
   [[nodiscard]] size_t size() const;

protected:
   VectorType& vector;
   const size_t start_index;
   const size_t end_index; // end_index is excluded
};

template <typename VectorType>
VectorView<VectorType>::VectorView(VectorType& vector, size_t start_index, size_t end_index):
      vector(vector), start_index(start_index), end_index(end_index) {
}

template <typename VectorType>
typename VectorType::value_type& VectorView<VectorType>::operator[](size_t index) {
   if (index < this->size()) {
      return this->vector[this->start_index + index];
   }
   throw std::runtime_error("You tried to access an element beyond the vector capacity.");
}

template <typename VectorType>
const typename VectorType::value_type& VectorView<VectorType>::operator[](size_t index) const {
   if (index < this->size()) {
      return this->vector[this->start_index + index];
   }
   throw std::runtime_error("You tried to access an element beyond the vector capacity.");
}

template <typename VectorType>
size_t VectorView<VectorType>::size() const {
   return this->end_index - this->start_index;
}

// free function
template <typename VectorType>
VectorView<VectorType> view(VectorType& vector, size_t start_index, size_t end_index) {
   return VectorView(vector, start_index, end_index);
}

#endif // UNO_VECTORVIEW_H