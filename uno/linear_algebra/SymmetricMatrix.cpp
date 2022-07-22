// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <exception>
#include <cassert>
#include <iomanip>
#include "SymmetricMatrix.hpp"

SymmetricMatrix::SymmetricMatrix(size_t dimension, size_t capacity) : dimension(dimension), capacity(capacity) {
   entries.reserve(capacity);
}

void SymmetricMatrix::reset() {
   this->number_nonzeros = 0;
}

double SymmetricMatrix::quadratic_product(const std::vector<double>& x, const std::vector<double>& y, size_t block_size) const {
   assert(x.size() == y.size() && "SymmetricMatrix::quadratic_product: the two vectors x and y do not have the same size");
   assert(block_size <= x.size() && "SymmetricMatrix::quadratic_product: the block size is larger than the vectors");

   double result = 0.;
   this->for_each([&](size_t i, size_t j, double entry) {
      if (i < block_size && j < block_size) {
         result += (i == j ? 1 : 2) * entry * x[i] * y[j];
      }
   });
   return result;
}

std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix& matrix) {
   stream << matrix.dimension << " variables, " << matrix.number_nonzeros << " non zeros:\n";
   matrix.print(stream);
   return stream;
}