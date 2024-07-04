// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICMATRIXFACTORY_H
#define UNO_SYMMETRICMATRIXFACTORY_H

#include "SymmetricMatrix.hpp"
#include "COOSymmetricMatrix.hpp"
#include "CSCSymmetricMatrix.hpp"

template <typename ElementType>
class SymmetricMatrixFactory {
public:
   static std::unique_ptr<SymmetricMatrix<ElementType>> create(const std::string& symmetric_matrix_type, size_t dimension, size_t capacity,
         bool use_regularization);
};

template <typename ElementType>
std::unique_ptr<SymmetricMatrix<ElementType>> SymmetricMatrixFactory<ElementType>::create(const std::string& symmetric_matrix_type, size_t dimension,
      size_t capacity, bool use_regularization) {
   if (symmetric_matrix_type == "COO") {
      return std::make_unique<COOSymmetricMatrix<ElementType>>(dimension, capacity, use_regularization);
   }
   else if (symmetric_matrix_type == "CSC") {
      return std::make_unique<CSCSymmetricMatrix<ElementType>>(dimension, capacity, use_regularization);
   }
   throw std::invalid_argument("Symmetric matrix type " + symmetric_matrix_type + " unknown");
}

#endif // UNO_SYMMETRICMATRIXFACTORY_H
