// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CSCSYMMETRICMATRIX_H
#define UNO_CSCSYMMETRICMATRIX_H

#include <ostream>
#include <cassert>
#include "SymmetricMatrix.hpp"
#include "tools/Infinity.hpp"

/*
 * Compressed Sparse Column
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 */

template <typename T>
class CSCSymmetricMatrix : public SymmetricMatrix<T> {
public:
   CSCSymmetricMatrix(size_t max_dimension, size_t original_capacity, bool use_regularization);

   void reset() override;
   void for_each(const std::function<void (size_t, size_t, T)>& f) const override;
   void for_each(size_t column_index, const std::function<void (size_t, T)>& f) const;
   void insert(T term, size_t row_index, size_t column_index) override;
   void finalize_column(size_t column_index) override;
   [[nodiscard]] T smallest_diagonal_entry() const override;
   void set_regularization(const std::function<T(size_t index)>& regularization_function) override;

   void print(std::ostream& stream) const override;

   static CSCSymmetricMatrix<T> identity(size_t dimension);

protected:
   // entries and row_indices have nnz elements
   // column_starts has dimension+1 elements
   std::vector<size_t> column_starts{};
   std::vector<size_t> row_indices{};
   size_t current_column{0};
   std::vector<T> diagonal_entries;
};

template <typename T>
CSCSymmetricMatrix<T>::CSCSymmetricMatrix(size_t max_dimension, size_t original_capacity, bool use_regularization):
      SymmetricMatrix<T>(max_dimension, original_capacity, use_regularization),
      column_starts(max_dimension + 1),
      diagonal_entries(max_dimension, T(0)) {
   this->entries.reserve(this->capacity);
   this->row_indices.reserve(this->capacity);
}

template <typename T>
void CSCSymmetricMatrix<T>::reset() {
   // empty the matrix
   SymmetricMatrix<T>::reset();
   this->row_indices.clear();
   initialize_vector<size_t>(this->column_starts, 0);
   this->current_column = 0;
   initialize_vector(this->diagonal_entries, T(0));
}

// generic iterator
template <typename T>
void CSCSymmetricMatrix<T>::for_each(const std::function<void (size_t, size_t, T)>& f) const {
   for (size_t j: Range(this->dimension)) {
      for (size_t k: Range(this->column_starts[j], this->column_starts[j + 1])) {
         const size_t i = this->row_indices[k];
         const T entry = this->entries[k];
         f(i, j, entry);
      }
   }
}

template <typename T>
void CSCSymmetricMatrix<T>::for_each(size_t column_index, const std::function<void (size_t, T)>& f) const {
   for (size_t k: Range(this->column_starts[column_index], this->column_starts[column_index + 1])) {
      const size_t i = this->row_indices[k];
      const T entry = this->entries[k];
      f(i, entry);
   }
}

template <typename T>
void CSCSymmetricMatrix<T>::insert(T term, size_t row_index, size_t column_index) {
   assert(column_index == this->current_column && "The previous columns should be finalized");

   this->entries.push_back(term);
   this->row_indices.push_back(row_index);
   this->column_starts[column_index + 1]++;
   this->number_nonzeros++;

   // possibly update diagonal
   if (row_index == column_index) {
      this->diagonal_entries[row_index] += term;
   }
}

template <typename T>
void CSCSymmetricMatrix<T>::finalize_column(size_t column_index) {
   assert(column_index == this->current_column && "You are not finalizing the current column");
   assert(column_index < this->dimension && "The dimension of the matrix was exceeded");

   // possibly add regularization
   if (this->use_regularization) {
      insert(0., column_index, column_index);
   }
   this->current_column++;

   // start the next column at the current start
   if (column_index < this->dimension - 1) {
      this->column_starts[column_index + 2] = this->column_starts[column_index + 1];
   }
}

template <typename T>
T CSCSymmetricMatrix<T>::smallest_diagonal_entry() const {
   T smallest_entry = INF<T>;
   for (size_t i: Range(this->dimension)) {
      smallest_entry = std::min(smallest_entry, this->diagonal_entries[i]);
   }
   return smallest_entry;
}

template <typename T>
void CSCSymmetricMatrix<T>::set_regularization(const std::function<T(size_t /*index*/)>& regularization_function) {
   assert(this->use_regularization && "You are trying to regularize a matrix where regularization was not preallocated.");

   for (size_t i: Range(this->dimension)) {
      // the regularization term is located at the end of the column, that is right before the start of the next column
      const size_t k = this->column_starts[i + 1] - 1;
      const T element = regularization_function(i);
      this->entries[k] = element;
      // update diagonal
      this->diagonal_entries[i] += element;
   }
}

template <typename T>
void CSCSymmetricMatrix<T>::print(std::ostream& stream) const {
   stream << "W = ";
   print_vector(stream, this->entries, 0, this->number_nonzeros);
   stream << "with column start: ";
   print_vector(stream, this->column_starts, 0, this->dimension + 1);
   stream << "and row index: ";
   print_vector(stream, this->row_indices, 0, this->number_nonzeros);
}

template <typename T>
CSCSymmetricMatrix<T> CSCSymmetricMatrix<T>::identity(size_t dimension) {
   CSCSymmetricMatrix<double> matrix(dimension, dimension, false);
   for (size_t i: Range(dimension)) {
      matrix.insert(1., i, i);
      matrix.finalize_column(i);
   }
   return matrix;
}

#endif // UNO_CSCSYMMETRICMATRIX_H