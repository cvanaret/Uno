#ifndef SYMMETRICMATRIX_H
#define SYMMETRICMATRIX_H

#include <vector>
#include <functional>
#include "SparseVector.hpp"

class SymmetricMatrix {
public:
   size_t dimension;
   size_t number_nonzeros{0};
   size_t capacity;

   SymmetricMatrix(size_t dimension, size_t capacity);
   virtual ~SymmetricMatrix() = default;

   void reset();
   [[nodiscard]] std::vector<double> product(const std::vector<double>& vector) const;
   [[nodiscard]] double quadratic_product(const std::vector<double>& x, const std::vector<double>& y) const;
   void add_outer_product(const SparseVector<double>& x, double scaling_factor = 1.);
   [[nodiscard]] double smallest_diagonal_entry() const;

   virtual void for_each(const std::function<void (size_t, size_t, double)>& f) const = 0;
   /* build the matrix incrementally */
   virtual void insert(double term, size_t row_index, size_t column_index) = 0;
   virtual void finalize(size_t column_index);
   virtual void add_identity_multiple(double multiple) = 0;
};

#endif // SYMMETRICMATRIX_H
