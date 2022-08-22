// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICMATRIXFACTORY_H
#define UNO_SYMMETRICMATRIXFACTORY_H

#include "SymmetricMatrix.hpp"
#include "COOSymmetricMatrix.hpp"
#include "CSCSymmetricMatrix.hpp"

template <typename T>
class SymmetricMatrixFactory {
public:
   static std::unique_ptr<SymmetricMatrix<T>> create(const std::string& symmetric_matrix_type, size_t dimension, size_t capacity,
         bool use_regularization);
};

template <typename T>
std::unique_ptr<SymmetricMatrix<T>> SymmetricMatrixFactory<T>::create(const std::string& symmetric_matrix_type, size_t dimension, size_t capacity,
      bool use_regularization) {
   if (symmetric_matrix_type == "COO") {
      return std::make_unique<COOSymmetricMatrix<T>>(dimension, capacity, use_regularization);
   }
   else if (symmetric_matrix_type == "CSC") {
      return std::make_unique<CSCSymmetricMatrix<T>>(dimension, capacity, use_regularization);
   }
   throw std::invalid_argument("Symmetric matrix type " + symmetric_matrix_type + " unknown");
}

#endif // UNO_SYMMETRICMATRIXFACTORY_H
