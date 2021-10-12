#ifndef COOSYMMETRICMATRIX_H
#define COOSYMMETRICMATRIX_H

#include "SymmetricMatrix.hpp"

class COOSymmetricMatrix : public SymmetricMatrix {
   /* Coordinate list */
public:
   std::vector<double> matrix;
   std::vector<size_t> row_indices;
   std::vector<size_t> column_indices;

   COOSymmetricMatrix(size_t dimension, size_t capacity);
   COOSymmetricMatrix(size_t dimension, std::vector<double> matrix, std::vector<size_t> row_indices, std::vector<size_t> column_indices);

   void for_each(const std::function<void (size_t, size_t, double)>& f) const override;
   void insert(double term, size_t row_index, size_t column_index) override;
   void add_identity_multiple(double multiple) override;

   friend std::ostream& operator<<(std::ostream& stream, const COOSymmetricMatrix& matrix);

protected:
   size_t find(size_t row_index, size_t column_index);
};

#endif // COOSYMMETRICMATRIX_H