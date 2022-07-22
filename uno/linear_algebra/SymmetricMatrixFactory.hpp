// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICMATRIXFACTORY_H
#define UNO_SYMMETRICMATRIXFACTORY_H

#include "SymmetricMatrix.hpp"

class SymmetricMatrixFactory {
public:
   static std::unique_ptr<SymmetricMatrix> create(const std::string& symmetric_matrix_type, size_t dimension, size_t capacity);
};

#endif // UNO_SYMMETRICMATRIXFACTORY_H
