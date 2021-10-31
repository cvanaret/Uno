#ifndef SYMMETRICMATRIX_H
#define SYMMETRICMATRIX_H

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
   [[nodiscard]] std::vector<double> product(const std::vector<double>& vector) const;
   [[nodiscard]] double quadratic_product(const std::vector<double>& x, const std::vector<double>& y) const;
   void add_outer_product(const SparseVector<double>& x, double column_entry = 1.);

   virtual void for_each(const std::function<void (size_t, size_t, double)>& f) const = 0;
   // build the matrix incrementally
   virtual void insert(double term, size_t row_index, size_t column_index) = 0;
   virtual void pop() = 0;
   virtual void finalize(size_t column_index);
   virtual void add_identity_multiple(double multiple) = 0;
   [[nodiscard]] virtual double smallest_diagonal_entry() const = 0;

   virtual void print(std::ostream& stream) const = 0;
   friend std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix& matrix);
};

#endif // SYMMETRICMATRIX_H
