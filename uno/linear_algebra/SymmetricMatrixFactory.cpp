// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#include <memory>
#include "SymmetricMatrixFactory.hpp"
#include "COOSymmetricMatrix.hpp"
#include "CSCSymmetricMatrix.hpp"

std::unique_ptr<SymmetricMatrix> SymmetricMatrixFactory::create(const std::string& symmetric_matrix_type, size_t dimension, size_t capacity) {
   if (symmetric_matrix_type == "COO") {
      return std::make_unique<COOSymmetricMatrix>(dimension, capacity);
   }
   else if (symmetric_matrix_type == "CSC") {
      const size_t padding_size = 1;
      return std::make_unique<CSCSymmetricMatrix>(dimension, capacity, padding_size);
   }
   throw std::invalid_argument("Symmetric matrix type unknown");
}