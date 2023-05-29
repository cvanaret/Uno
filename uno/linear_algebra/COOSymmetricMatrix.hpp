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
template <typename T>
class COOSymmetricMatrix: public SymmetricMatrix<T> {
   // Coordinate list
public:
   COOSymmetricMatrix(size_t max_dimension, size_t original_capacity, bool use_regularization);

   void reset() override;
   void for_each(const std::function<void (size_t, size_t, T)>& f) const override;
   void insert(T term, size_t row_index, size_t column_index) override;
   void finalize_column(size_t column_index) override;
   [[nodiscard]] T smallest_diagonal_entry() const override;
   void set_regularization(const std::function<T(size_t index)>& regularization_function) override;

   void print(std::ostream& stream) const override;

   static COOSymmetricMatrix<T> zero(size_t dimension);

protected:
   std::vector<size_t> row_indices;
   std::vector<size_t> column_indices;
   std::vector<T> diagonal_entries;

   void initialize_regularization();
};

// implementation

template <typename T>
COOSymmetricMatrix<T>::COOSymmetricMatrix(size_t max_dimension, size_t original_capacity, bool use_regularization):
      SymmetricMatrix<T>(max_dimension, original_capacity, use_regularization),
      diagonal_entries(max_dimension, T(0)) {
   this->row_indices.reserve(this->capacity);
   this->column_indices.reserve(this->capacity);

   // initialize regularization terms
   if (this->use_regularization) {
      this->initialize_regularization();
   }
}

template <typename T>
void COOSymmetricMatrix<T>::reset() {
   // empty the matrix
   SymmetricMatrix<T>::reset();
   this->row_indices.clear();
   this->column_indices.clear();
   initialize_vector(this->diagonal_entries, T(0));

   // initialize regularization terms
   if (this->use_regularization) {
      this->initialize_regularization();
   }
}

// generic iterator
template <typename T>
void COOSymmetricMatrix<T>::for_each(const std::function<void(size_t, size_t, T)>& f) const {
   for (size_t k: Range(this->number_nonzeros)) {
      const size_t i = this->row_indices[k];
      const size_t j = this->column_indices[k];
      const T entry = this->entries[k];
      f(i, j, entry);
   }
}

template <typename T>
void COOSymmetricMatrix<T>::insert(T term, size_t row_index, size_t column_index) {
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

template <typename T>
void COOSymmetricMatrix<T>::finalize_column(size_t /*column_index*/) {
   // do nothing
}

template <typename T>
T COOSymmetricMatrix<T>::smallest_diagonal_entry() const {
   T smallest_entry = INF<T>;
   for (size_t i: Range(this->dimension)) {
      smallest_entry = std::min(smallest_entry, this->diagonal_entries[i]);
   }
   return smallest_entry;
}

template <typename T>
void COOSymmetricMatrix<T>::set_regularization(const std::function<T(size_t /*index*/)>& regularization_function) {
   assert(this->use_regularization && "You are trying to regularize a matrix where regularization was not preallocated.");

   // the regularization terms (that lie at the start of the entries vector) can be directly modified
   for (size_t i: Range(this->dimension)) {
      const T element = regularization_function(i);
      this->entries[i] = element;
      // update diagonal
      this->diagonal_entries[i] += element;
   }
}

template <typename T>
void COOSymmetricMatrix<T>::print(std::ostream& stream) const {
   this->for_each([&](size_t i, size_t j, T entry) {
      stream << "m(" << i << ", " << j << ") = " << entry << '\n';
   });
}

template <typename T>
void COOSymmetricMatrix<T>::initialize_regularization() {
   // introduce elements at the start of the entries
   for (size_t i: Range(this->dimension)) {
      this->insert(T(0), i, i);
   }
}

template <typename T>
COOSymmetricMatrix<T> COOSymmetricMatrix<T>::zero(size_t dimension) {
   return COOSymmetricMatrix<T>(dimension, 0, false);
}

#endif // UNO_COOSYMMETRICMATRIX_H