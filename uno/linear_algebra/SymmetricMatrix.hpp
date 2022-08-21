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
   void product(const std::vector<T>& vector, std::vector<T>& result) const;
   [[nodiscard]] T quadratic_product(const std::vector<T>& x, const std::vector<T>& y, size_t block_size) const;

   virtual void for_each(const std::function<void (size_t, size_t, T)>& f) const = 0;
   // build the matrix incrementally
   virtual void insert(T term, size_t row_index, size_t column_index) = 0;
   // this method will be used by the CSCSymmetricMatrix subclass
   virtual void finalize_column(size_t column_index) = 0;
   [[nodiscard]] virtual T smallest_diagonal_entry() const = 0;
   virtual void set_regularization(const std::function<T(size_t index)>& regularization_function) = 0;
   [[nodiscard]] const T* raw_pointer() const;

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
void SymmetricMatrix<T>::product(const std::vector<T>& vector, std::vector<T>& result) const {
   assert(this->dimension == vector.size() && "The matrix and the vector do not have the same size");

   initialize_vector(result, T(0));
   this->for_each([&](size_t i, size_t j, T entry) {
      result[i] += entry * vector[j];
      // off-diagonal term
      if (i != j) {
         result[j] += entry * vector[i];
      }
   });
}

template <typename T>
T SymmetricMatrix<T>::quadratic_product(const std::vector<T>& x, const std::vector<T>& y, size_t block_size) const {
   assert(x.size() == y.size() && "SymmetricMatrix::quadratic_product: the two vectors x and y do not have the same size");
   assert(block_size <= x.size() && "SymmetricMatrix::quadratic_product: the block size is larger than the vectors");

   T result = T(0);
   this->for_each([&](size_t i, size_t j, T entry) {
      if (i < block_size && j < block_size) {
         result += (i == j ? T(1) : T(2)) * entry * x[i] * y[j];
      }
   });
   return result;
}

template <typename T>
const T* SymmetricMatrix<T>::raw_pointer() const {
   return this->entries.data();
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix<T>& matrix) {
   stream << "Dimension: " << matrix.dimension << ", number of nonzeros: " << matrix.number_nonzeros << '\n';
   matrix.print(stream);
   return stream;
}

#endif // UNO_SYMMETRICMATRIX_H
