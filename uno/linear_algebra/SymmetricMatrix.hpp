// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICMATRIX_H
#define UNO_SYMMETRICMATRIX_H

#include <vector>
#include <functional>
#include <cassert>
#include "Matrix.hpp"
#include "Vector.hpp"
#include "SparseVector.hpp"

// abstract class
template <typename ElementType>
class SymmetricMatrix: public Matrix<ElementType> {
public:
   using value_type = ElementType;

   size_t dimension;
   size_t number_nonzeros{0};
   size_t capacity;

   SymmetricMatrix(size_t max_dimension, size_t original_capacity, bool use_regularization);
   virtual ~SymmetricMatrix() = default;

   virtual void reset();

   [[nodiscard]] size_t number_rows() const;
   [[nodiscard]] size_t number_columns() const;

   template <typename VectorType, typename ResultType>
   void product(const VectorType& vector, ResultType& result) const;

   ElementType quadratic_product(const std::vector<ElementType>& x, const std::vector<ElementType>& y) const;

   // build the matrix incrementally
   virtual void insert(ElementType term, size_t row_index, size_t column_index) = 0;
   // this method will be used by the CSCSymmetricMatrix subclass
   virtual void finalize_column(size_t column_index) = 0;
   [[nodiscard]] virtual ElementType smallest_diagonal_entry() const = 0;
   virtual void set_regularization(const std::function<ElementType(size_t /*index*/)>& regularization_function) = 0;

   [[nodiscard]] const ElementType* data_raw_pointer() const;

   virtual void print(std::ostream& stream) const = 0;
   template <typename U>
   friend std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix<U>& matrix);

protected:
   std::vector<ElementType> entries{};
   // regularization
   const bool use_regularization;
};

// implementation

template <typename ElementType>
SymmetricMatrix<ElementType>::SymmetricMatrix(size_t max_dimension, size_t original_capacity, bool use_regularization):
      Matrix<ElementType>(),
      dimension(max_dimension),
      // if regularization is used, allocate the necessary space
      capacity(original_capacity + (use_regularization ? max_dimension : 0)),
      use_regularization(use_regularization) {
   this->entries.reserve(this->capacity);
}

template <typename ElementType>
void SymmetricMatrix<ElementType>::reset() {
   this->number_nonzeros = 0;
   this->entries.clear();
}

template <typename ElementType>
size_t SymmetricMatrix<ElementType>::number_rows() const {
   return this->dimension;
}

template <typename ElementType>
size_t SymmetricMatrix<ElementType>::number_columns() const {
   return this->dimension;
}

template <typename ElementType>
template <typename VectorType, typename ResultType>
void SymmetricMatrix<ElementType>::product(const VectorType& vector, ResultType& result) const {
   this->for_each([&](size_t row_index, size_t column_index, double entry) {
      result[row_index] += entry*vector[column_index];
      // off-diagonal term on the other side of the diagonal
      if (row_index != column_index) {
         result[column_index] += entry*vector[row_index];
      }
   });
}

template <typename ElementType>
ElementType SymmetricMatrix<ElementType>::quadratic_product(const std::vector<ElementType>& x, const std::vector<ElementType>& y) const {
   assert(x.size() == y.size() && "SymmetricMatrix::quadratic_product: the two vectors x and y do not have the same size");

   ElementType result = ElementType(0);
   this->for_each([&](size_t row_index, size_t column_index, ElementType entry) {
      result += (row_index == column_index ? ElementType(1) : ElementType(2)) * entry * x[row_index] * y[column_index];
   });
   return result;
}

template <typename ElementType>
const ElementType* SymmetricMatrix<ElementType>::data_raw_pointer() const {
   return this->entries.data();
}

template <typename ElementType>
std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix<ElementType>& matrix) {
   stream << "Dimension: " << matrix.dimension << ", number of nonzeros: " << matrix.number_nonzeros << '\n';
   matrix.print(stream);
   return stream;
}

#endif // UNO_SYMMETRICMATRIX_H
