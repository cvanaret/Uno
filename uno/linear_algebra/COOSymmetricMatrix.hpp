// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOSYMMETRICMATRIX_H
#define UNO_COOSYMMETRICMATRIX_H

#include <ostream>
#include <cassert>
#include "SymmetricMatrix.hpp"
#include "tools/Infinity.hpp"

/*
 * Coordinate list
 * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
 */
template <typename ElementType>
class COOSymmetricMatrix: public SymmetricMatrix<ElementType> {
   // Coordinate list
public:
   COOSymmetricMatrix(size_t max_dimension, size_t original_capacity, bool use_regularization);

   void reset() override;
   void for_each(const std::function<void(size_t, size_t, ElementType)>& f) const override;
   void insert(ElementType term, size_t row_index, size_t column_index) override;
   void finalize_column(size_t column_index) override;
   [[nodiscard]] ElementType smallest_diagonal_entry() const override;
   void set_regularization(const std::function<ElementType(size_t index)>& regularization_function) override;

   void print(std::ostream& row_index) const override;

   static COOSymmetricMatrix<ElementType> zero(size_t dimension);

protected:
   std::vector<size_t> row_indices;
   std::vector<size_t> column_indices;
   std::vector<ElementType> diagonal_entries;

   void initialize_regularization();
};

// implementation

template <typename ElementType>
COOSymmetricMatrix<ElementType>::COOSymmetricMatrix(size_t max_dimension, size_t original_capacity, bool use_regularization):
      SymmetricMatrix<ElementType>(max_dimension, original_capacity, use_regularization),
      diagonal_entries(max_dimension, ElementType(0)) {
   this->row_indices.reserve(this->capacity);
   this->column_indices.reserve(this->capacity);

   // initialize regularization terms
   if (this->use_regularization) {
      this->initialize_regularization();
   }
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::reset() {
   // empty the matrix
   SymmetricMatrix<ElementType>::reset();
   this->row_indices.clear();
   this->column_indices.clear();
   initialize_vector(this->diagonal_entries, ElementType(0));

   // initialize regularization terms
   if (this->use_regularization) {
      this->initialize_regularization();
   }
}

// generic iterator
template <typename ElementType>
void COOSymmetricMatrix<ElementType>::for_each(const std::function<void(size_t, size_t, ElementType)>& f) const {
   for (size_t k: Range(this->number_nonzeros)) {
      const size_t row_index = this->row_indices[k];
      const size_t column_index = this->column_indices[k];
      const ElementType entry = this->entries[k];
      f(row_index, column_index, entry);
   }
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::insert(ElementType term, size_t row_index, size_t column_index) {
   assert(this->number_nonzeros <= this->row_indices.size() && "The COO matrix doesn't have a sufficient capacity");
   this->entries.push_back(term);
   this->row_indices.push_back(row_index);
   this->column_indices.push_back(column_index);
   this->number_nonzeros++;

   // possibly update diagonal
   if (row_index == column_index) {
      this->diagonal_entries[row_index] += term;
   }
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::finalize_column(size_t /*column_index*/) {
   // do nothing
}

template <typename ElementType>
ElementType COOSymmetricMatrix<ElementType>::smallest_diagonal_entry() const {
   ElementType smallest_entry = INF<ElementType>;
   for (size_t row_index: Range(this->dimension)) {
      smallest_entry = std::min(smallest_entry, this->diagonal_entries[row_index]);
   }
   return smallest_entry;
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::set_regularization(const std::function<ElementType(size_t /*index*/)>& regularization_function) {
   assert(this->use_regularization && "You are trying to regularize a matrix where regularization was not preallocated.");

   // the regularization terms (that lie at the start of the entries vector) can be directly modified
   for (size_t row_index: Range(this->dimension)) {
      const ElementType element = regularization_function(row_index);
      this->entries[row_index] = element;
      // update diagonal
      this->diagonal_entries[row_index] += element;
   }
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::print(std::ostream& stream) const {
   this->for_each([&](size_t row_index, size_t column_index, ElementType entry) {
      stream << "m(" << row_index << ", " << column_index << ") = " << entry << '\n';
   });
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::initialize_regularization() {
   // introduce elements at the start of the entries
   for (size_t row_index: Range(this->dimension)) {
      this->insert(ElementType(0), row_index, row_index);
   }
}

template <typename ElementType>
COOSymmetricMatrix<ElementType> COOSymmetricMatrix<ElementType>::zero(size_t dimension) {
   return COOSymmetricMatrix<ElementType>(dimension, 0, false);
}

#endif // UNO_COOSYMMETRICMATRIX_H