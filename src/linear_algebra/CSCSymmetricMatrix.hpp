#ifndef CSCSYMMETRICMATRIX_H
#define CSCSYMMETRICMATRIX_H

#include "SymmetricMatrix.hpp"
#include "COOSymmetricMatrix.hpp"

class CSCSymmetricMatrix : public SymmetricMatrix {
   /* Compressed Sparse Column */
public:
   std::vector<double> matrix;
   std::vector<int> column_start;
   std::vector<int> row_index;

   CSCSymmetricMatrix(int dimension, size_t maximum_number_nonzeros);
   CSCSymmetricMatrix(std::vector<double> matrix, const std::vector<int>& column_start, std::vector<int> row_number, int capacity);

   void for_each(const std::function<void (int, int, double)>& f) const override;
   void insert(double term, int row_index, int column_index) override;
   static CSCSymmetricMatrix identity(int dimension);

   void add_identity_multiple(double multiple) override;
   COOSymmetricMatrix to_COO();

   friend std::ostream& operator<<(std::ostream& stream, CSCSymmetricMatrix& matrix);
   friend std::ostream& operator<<(std::ostream& stream, const CSCSymmetricMatrix& matrix);
};

#endif // CSCSYMMETRICMATRIX_H