#ifndef CSCSYMMETRICMATRIX_H
#define CSCSYMMETRICMATRIX_H

#include "SymmetricMatrix.hpp"
#include "COOSymmetricMatrix.hpp"

enum CSCDiagonal {
   SPARSE,
   EXPLICIT
};

class CSCSymmetricMatrix : public SymmetricMatrix {
   /* Compressed Sparse Column */
public:
   std::vector<double> matrix;
   std::vector<int> column_start;
   std::vector<int> row_index;

   CSCSymmetricMatrix(size_t dimension, size_t maximum_number_nonzeros, size_t padding_size = 0);
   CSCSymmetricMatrix(std::vector<double> matrix, const std::vector<int>& column_start, std::vector<int> row_number, int capacity);

   void for_each(const std::function<void (int, int, double)>& f) const override;
   void for_each(size_t column_index, const std::function<void (int, double)>& f) const;
   void insert(double term, size_t row_index, size_t column_index) override;
   void finalize(size_t column_index);
   void add_identity_multiple(double factor) override;
   void force_explicit_diagonal_elements();
   void remove_variables(const std::vector<int>& variable_indices);
   COOSymmetricMatrix to_COO();

   static CSCSymmetricMatrix identity(int dimension);

   friend std::ostream& operator<<(std::ostream& stream, CSCSymmetricMatrix& matrix);
   friend std::ostream& operator<<(std::ostream& stream, const CSCSymmetricMatrix& matrix);

protected:
   // when elements are inserted one by one, keep track of the current column
   size_t current_column{0};
   std::vector<size_t> remaining_column_padding;
};

#endif // CSCSYMMETRICMATRIX_H