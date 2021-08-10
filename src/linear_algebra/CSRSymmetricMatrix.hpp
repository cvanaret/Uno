#ifndef CSRSYMMETRICMATRIX_H
#define CSRSYMMETRICMATRIX_H

#include "SymmetricMatrix.hpp"
#include "COOSymmetricMatrix.hpp"

class CSRSymmetricMatrix : public SymmetricMatrix {
   /* Compressed Sparse Row */
public:
   std::vector<double> matrix;
   std::vector<size_t> row_start;
   std::vector<size_t> column_index;

   CSRSymmetricMatrix(size_t dimension, size_t maximum_number_nonzeros);
   CSRSymmetricMatrix(const std::vector<double>& matrix, const std::vector<size_t>& column_start, const std::vector<size_t>& row_number, size_t capacity);

   void for_each(const std::function<void (size_t, size_t, double)>& f) const override;
   void insert(double term, size_t row_index, size_t column_index) override;
   static CSRSymmetricMatrix identity(size_t dimension);

   CSRSymmetricMatrix add_identity_multiple(double multiple);
   COOSymmetricMatrix to_COO();

   friend std::ostream& operator<<(std::ostream& stream, CSRSymmetricMatrix& matrix);
   friend std::ostream& operator<<(std::ostream& stream, const CSRSymmetricMatrix& matrix);
};

#endif // CSRSYMMETRICMATRIX_H