#ifndef COOSYMMETRICMATRIX_H
#define COOSYMMETRICMATRIX_H

#include "SymmetricMatrix.hpp"

class COOSymmetricMatrix : public SymmetricMatrix {
   /* Coordinate list */
public:
   std::vector<double> matrix;
   std::vector<int> row_indices;
   std::vector<int> column_indices;

   COOSymmetricMatrix(size_t dimension, size_t capacity);
   COOSymmetricMatrix(size_t dimension, std::vector<double> matrix, std::vector<int> row_indices, std::vector<int> column_indices);

   void for_each(const std::function<void (int, int, double)>& f) const override;
   void insert(double term, size_t row_index, size_t column_index) override;
   void add_identity_multiple(double multiple) override;

   friend std::ostream& operator<<(std::ostream& stream, COOSymmetricMatrix& matrix);
   friend std::ostream& operator<<(std::ostream& stream, const COOSymmetricMatrix& matrix);
};

#endif // COOSYMMETRICMATRIX_H