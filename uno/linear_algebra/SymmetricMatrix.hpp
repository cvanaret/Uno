// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#ifndef UNO_SYMMETRICMATRIX_H
#define UNO_SYMMETRICMATRIX_H

#include <vector>
#include <functional>
#include "SparseVector.hpp"

class SymmetricMatrix {
public:
   std::vector<double> entries;
   size_t dimension;
   size_t number_nonzeros{0};
   size_t capacity;

   SymmetricMatrix(size_t dimension, size_t capacity);
   virtual ~SymmetricMatrix() = default;

   virtual void reset();
   [[nodiscard]] double quadratic_product(const std::vector<double>& x, const std::vector<double>& y, size_t block_size) const;

   virtual void for_each(const std::function<void (size_t, size_t, double)>& f) const = 0;
   // build the matrix incrementally
   virtual void insert(double term, size_t row_index, size_t column_index) = 0;
   virtual void pop() = 0;
   virtual void finalize(size_t column_index) = 0;
   virtual void add_identity_multiple(double multiple, size_t number_variables) = 0;
   [[nodiscard]] virtual double smallest_diagonal_entry() const = 0;

   virtual void print(std::ostream& stream) const = 0;
   friend std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix& matrix);
};

#endif // UNO_SYMMETRICMATRIX_H
