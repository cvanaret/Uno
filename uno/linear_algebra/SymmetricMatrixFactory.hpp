#ifndef SYMMETRICMATRIXFACTORY_H
#define SYMMETRICMATRIXFACTORY_H

#include "SymmetricMatrix.hpp"

class SymmetricMatrixFactory {
public:
   static std::unique_ptr<SymmetricMatrix> create(const std::string& symmetric_matrix_type, size_t dimension, size_t capacity);
};

#endif // SYMMETRICMATRIXFACTORY_H
