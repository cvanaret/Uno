// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICMATRIX_H
#define UNO_SYMMETRICMATRIX_H

#include <vector>
#include <functional>
#include <cassert>
#include "Vector.hpp"
#include "SparseVector.hpp"

template <typename T>
class SymmetricMatrix {
public:
   size_t dimension;
   size_t number_nonzeros{0};
   size_t capacity;

   SymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization);
   virtual ~SymmetricMatrix() = default;

   virtual void reset();

   virtual void for_each(const std::function<void (size_t, size_t, T)>& f) const = 0;
   // build the matrix incrementally
   virtual void insert(T term, size_t row_index, size_t column_index) = 0;
   // this method will be used by the CSCSymmetricMatrix subclass
   virtual void finalize_column(size_t column_index) = 0;
   [[nodiscard]] virtual T smallest_diagonal_entry() const = 0;
   virtual void set_regularization(const std::function<T(size_t index)>& regularization_function) = 0;
   [[nodiscard]] const T* data_raw_pointer() const;

   virtual void print(std::ostream& stream) const = 0;
   template <typename U>
   friend std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix<U>& matrix);

protected:
   std::vector<T> entries{};
   // regularization
   const bool use_regularization;
};

// implementation

template <typename T>
SymmetricMatrix<T>::SymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization) :
      dimension(dimension),
      // if regularization is used, allocate the necessary space
      capacity(original_capacity + (use_regularization ? dimension : 0)),
      use_regularization(use_regularization) {
   this->entries.reserve(this->capacity);
}

template <typename T>
void SymmetricMatrix<T>::reset() {
   this->number_nonzeros = 0;
   this->entries.clear();
}

template <typename T>
const T* SymmetricMatrix<T>::data_raw_pointer() const {
   return this->entries.data();
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix<T>& matrix) {
   stream << "Dimension: " << matrix.dimension << ", number of nonzeros: " << matrix.number_nonzeros << '\n';
   matrix.print(stream);
   return stream;
}

#endif // UNO_SYMMETRICMATRIX_H
