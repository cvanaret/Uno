// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MATRIX_H
#define UNO_MATRIX_H

namespace uno {
   // abstract class
   template <typename IndexType, typename ElementType>
   class Matrix {
   public:
      Matrix(size_t number_rows, size_t number_columns, size_t number_nonzeros);
      virtual ~Matrix() = default;

      [[nodiscard]] virtual size_t get_number_rows() const;
      [[nodiscard]] virtual size_t get_number_columns() const;
      [[nodiscard]] virtual size_t get_number_nonzeros() const;

      virtual void insert(ElementType term, IndexType row_index, IndexType column_index) = 0;
      virtual void finalize_column(IndexType column_index) = 0;

   protected:
      const size_t number_rows, number_columns, number_nonzeros;
   };

   template <typename IndexType, typename ElementType>
   Matrix<IndexType, ElementType>::Matrix(size_t number_rows, size_t number_columns, size_t number_nonzeros):
      number_rows(number_rows), number_columns(number_columns), number_nonzeros(number_nonzeros) { }

   template <typename IndexType, typename ElementType>
   size_t Matrix<IndexType, ElementType>::get_number_rows() const {
      return this->number_rows;
   }

   template <typename IndexType, typename ElementType>
   size_t Matrix<IndexType, ElementType>::get_number_columns() const {
      return this->number_columns;
   }

   template <typename IndexType, typename ElementType>
   size_t Matrix<IndexType, ElementType>::get_number_nonzeros() const {
      return this->number_nonzeros;
   }
} // namespace

#endif // UNO_MATRIX_H