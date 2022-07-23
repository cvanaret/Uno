// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <memory>
#include "SymmetricMatrixFactory.hpp"
#include "COOSymmetricMatrix.hpp"
#include "CSCSymmetricMatrix.hpp"

std::unique_ptr<SymmetricMatrix> SymmetricMatrixFactory::create(const std::string& symmetric_matrix_type, size_t dimension, size_t capacity,
      bool use_regularization) {
   if (symmetric_matrix_type == "COO") {
      return std::make_unique<COOSymmetricMatrix>(dimension, capacity, use_regularization);
   }
   else if (symmetric_matrix_type == "CSC") {
      return std::make_unique<CSCSymmetricMatrix>(dimension, capacity, use_regularization);
   }
   throw std::invalid_argument("Symmetric matrix type unknown");
}