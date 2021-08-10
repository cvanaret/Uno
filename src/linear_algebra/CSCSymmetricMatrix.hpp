#ifndef CSCSYMMETRICMATRIX_H
#define CSCSYMMETRICMATRIX_H

#include "SymmetricMatrix.hpp"
#include "COOSymmetricMatrix.hpp"

class CSCSymmetricMatrix : public SymmetricMatrix {
   /* Compressed Sparse Column */
public:
   std::vector<double> matrix;
   std::vector<size_t> column_start;
   std::vector<size_t> row_index;

   CSCSymmetricMatrix(size_t dimension, size_t maximum_number_nonzeros);
   CSCSymmetricMatrix(const std::vector<double>& matrix, const std::vector<size_t>& column_start, const std::vector<size_t>& row_number, size_t capacity);

   void for_each(const std::function<void (size_t, size_t, double)>& f) const override;
   void insert(double term, size_t row_index, size_t column_index) override;
   static CSCSymmetricMatrix identity(size_t dimension);

   CSCSymmetricMatrix add_identity_multiple(double multiple);
   COOSymmetricMatrix to_COO();

   friend std::ostream& operator<<(std::ostream& stream, CSCSymmetricMatrix& matrix);
   friend std::ostream& operator<<(std::ostream& stream, const CSCSymmetricMatrix& matrix);
};

#endif // CSCSYMMETRICMATRIX_H