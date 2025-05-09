// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PARAMETERIZEDMATRIX_H
#define UNO_PARAMETERIZEDMATRIX_H

#include <stdexcept>
#include "linear_algebra/COOSparseStorage.hpp"
#include "linear_algebra/CSCSparseStorage.hpp"

enum class Storage {COO, CSC};

enum class MatrixType {RECTANGULAR, SQUARE, SYMMETRIC};

namespace uno {
   template <Storage sparse_storage, MatrixType matrix_type>
   class SparseMatrix {
   public:
      SparseMatrix(size_t number_rows, size_t number_columns):
            // bullshit initialization
            storage(number_rows, 0, false) {
         if constexpr (matrix_type == MatrixType::SQUARE || matrix_type == MatrixType::SYMMETRIC) {
            if (number_rows != number_columns) {
               throw std::runtime_error("The square matrix has different dimensions");
            }
         }
      }

   protected:
      // set the type of the sparse storage at compile time
      std::conditional_t<sparse_storage == Storage::COO,
         COOSparseStorage<size_t, double>,
         CSCSparseStorage<size_t, double>> storage;
   };
} // namespace

#endif // UNO_PARAMETERIZEDMATRIX_H