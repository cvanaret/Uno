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

   [[nodiscard]] std::vector<double> product(const std::vector<double>& vector) const;
   [[nodiscard]] double quadratic_product(const std::vector<double>& x, const std::vector<double>& y) const;
   void add_outer_product(const SparseVector& x, double scaling_factor = 1.);
   [[nodiscard]] double smallest_diagonal_entry() const;

   virtual void for_each(const std::function<void (size_t, size_t, double)>& f) const = 0;
   /* build the matrix incrementally */
   virtual void insert(double term, size_t row_index, size_t column_index) = 0;
};

//class UnoMatrix : public SymmetricMatrix {
//   /* Coordinate list */
//public:
//   UnoMatrix(size_t dimension, size_t number_nonzeros);
//
//   SparseVector matrix;
//
//   void insert(double term, size_t row_index, size_t column_index) override;
//   std::vector<double> product(const std::vector<double>& vector) override;
//   void add_matrix(UnoMatrix& other_matrix, double factor);
//
//   double norm_1();
//   COOSymmetricMatrix to_COO();
//   COOSymmetricMatrix to_COO(const std::unordered_map<size_t, size_t>& mask);
//   //CSCSymmetricMatrix to_CSC();
//   //CSCSymmetricMatrix to_CSC(const std::unordered_map<int, int>& mask);
//
//   friend std::ostream& operator<<(std::ostream& stream, UnoMatrix& matrix);
//};

#endif // SYMMETRICMATRIX_H
