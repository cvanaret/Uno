// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICMATRIX_H
#define UNO_SYMMETRICMATRIX_H

#include <vector>
#include <functional>
#include "SparseVector.hpp"

class SymmetricMatrix {
public:
   size_t dimension;
   size_t number_nonzeros{0};
   size_t capacity;

   SymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization);
   virtual ~SymmetricMatrix() = default;

   virtual void reset();
   [[nodiscard]] double quadratic_product(const std::vector<double>& x, const std::vector<double>& y, size_t block_size) const;

   virtual void for_each(const std::function<void (size_t, size_t, double)>& f) const = 0;
   // build the matrix incrementally
   virtual void insert(double term, size_t row_index, size_t column_index) = 0;
   virtual void finalize(size_t column_index) = 0;
   [[nodiscard]] virtual double smallest_diagonal_entry() const = 0;
   virtual void set_regularization(const std::function<double(size_t index)>& f) = 0;
   [[nodiscard]] const double* raw_pointer() const;

   virtual void print(std::ostream& stream) const = 0;
   friend std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix& matrix);

protected:
   std::vector<double> entries{};
   // regularization
   const bool use_regularization;
};

#endif // UNO_SYMMETRICMATRIX_H
