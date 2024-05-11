// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RECTANGULARMATRIX_H
#define UNO_RECTANGULARMATRIX_H

#include "SparseVector.hpp"
#include "Matrix.hpp"

// TODO use more appropriate sparse representation

template <typename ElementType>
using RectangularMatrix = std::vector<SparseVector<ElementType>>;

/*
// row major matrix
template <typename ElementType>
class RectangularMatrix: public Matrix<ElementType> {
public:
   class iterator {
   public:
      iterator(const RectangularMatrix<ElementType>& matrix, size_t row_index): matrix(matrix), row_index(row_index) { }
      const SparseVector<ElementType>& operator*() const { return this->matrix[this->row_index]; }
      // prefix increment
      iterator& operator++() {
         this->row_index++;
         return *this;
      }

      friend bool operator==(const iterator& a, const iterator& b) { return &a.matrix == &b.matrix && a.row_index == b.row_index; };
      friend bool operator!=(const iterator& a, const iterator& b) { return &a.matrix != &b.matrix || a.row_index != b.row_index; };

   protected:
      const RectangularMatrix<ElementType>& matrix;
      size_t row_index;
   };

   RectangularMatrix(size_t number_rows, size_t number_columns):
         matrix(number_rows), number_of_rows(number_rows), number_of_columns(number_columns) {
      for (auto& constraint_gradient: this->matrix) {
         constraint_gradient.reserve(number_columns);
      }
   }

   RectangularMatrix(RectangularMatrix<ElementType>&& other_matrix) noexcept:
      matrix(std::move(other_matrix.matrix)), number_of_rows(other_matrix.number_of_rows), number_of_columns(other_matrix.number_of_columns) { }

   RectangularMatrix& operator=(RectangularMatrix<ElementType>&& other_matrix)  noexcept {
      if (this != &other_matrix) {
         this->matrix = std::move(other_matrix.matrix);
      }
      return *this;
   }

   [[nodiscard]] size_t number_rows() const { return this->number_of_rows; }
   [[nodiscard]] size_t number_columns() const { return this->number_of_columns; }

   SparseVector<ElementType>& operator[](size_t row_index) { return this->matrix[row_index]; }
   const SparseVector<ElementType>& operator[](size_t row_index) const { return this->matrix[row_index]; }

   void insert(size_t row_index, size_t column_index, ElementType element) {
      this->matrix[row_index].insert(column_index, element);
   }

   void for_each(const std::function<void(size_t, size_t, ElementType)>& f) const { }

   iterator begin() const { return iterator(*this, 0); }
   iterator end() const { return iterator(*this, this->number_rows()); }

protected:
   std::vector<SparseVector<ElementType>> matrix;
   const size_t number_of_rows;
   const size_t number_of_columns;
};
*/

#endif // UNO_RECTANGULARMATRIX_H