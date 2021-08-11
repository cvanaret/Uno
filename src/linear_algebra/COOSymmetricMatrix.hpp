#ifndef COOSYMMETRICMATRIX_H
#define COOSYMMETRICMATRIX_H

#include "SymmetricMatrix.hpp"

class COOSymmetricMatrix : public SymmetricMatrix {
   /* Coordinate list */
public:
   std::vector<double> matrix;
   std::vector<int> row_indices;
   std::vector<int> column_indices;

   COOSymmetricMatrix(int dimension, size_t capacity);

   void for_each(const std::function<void (int, int, double)>& f) const override;
   void insert(double term, int row_index, int column_index) override;

   double norm_1();
   /*COOSymmetricMatrix add_identity_multiple(double multiple);*/

   friend std::ostream& operator<<(std::ostream& stream, COOSymmetricMatrix& matrix);
   friend std::ostream& operator<<(std::ostream& stream, const COOSymmetricMatrix& matrix);
};

#endif // COOSYMMETRICMATRIX_H