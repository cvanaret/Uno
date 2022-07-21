// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#include <memory>
#include "SymmetricMatrixFactory.hpp"
#include "COOSymmetricMatrix.hpp"

std::unique_ptr<SymmetricMatrix> SymmetricMatrixFactory::create(const std::string& symmetric_matrix_type, size_t dimension, size_t capacity) {
   if (symmetric_matrix_type == "COO") {
      return std::make_unique<COOSymmetricMatrix>(dimension, capacity);
   }
   throw std::invalid_argument("Symmetric matrix type unknown");
}