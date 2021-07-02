#ifndef MATRIX_H
#define MATRIX_H

#include <ostream>
#include <map>
#include <vector>
#include "SparseVector.hpp"

class Matrix {
public:
   Matrix(size_t dimension, size_t capacity, short fortran_indexing);
   virtual ~Matrix() = default;

   size_t dimension;
   size_t number_nonzeros;
   size_t capacity;
   short fortran_indexing;

   /* build the matrix incrementally */
   virtual void insert(double term, size_t row_index, size_t column_index) = 0;
   virtual std::vector<double> product(const std::vector<double>& vector) = 0;

   virtual double quadratic_product(const std::vector<double>& x, const std::vector<double>& y);
   void add_outer_product(const SparseVector& x, double scaling_factor = 1.);
};

class COOMatrix : public Matrix {
   /* Coordinate list */
public:
   COOMatrix(size_t dimension, size_t capacity, short fortran_indexing);

   std::vector<double> matrix;
   std::vector<int> row_indices;
   std::vector<int> column_indices;

   void insert(double term, size_t row_index, size_t column_index) override;
   std::vector<double> product(const std::vector<double>& vector) override;

   double norm_1();
   /*COOMatrix add_identity_multiple(double multiple);*/

   friend std::ostream& operator<<(std::ostream& stream, COOMatrix& matrix);
   friend std::ostream& operator<<(std::ostream& stream, const COOMatrix& matrix);
};

// forward declaration
class UnoMatrix;

class CSCMatrix : public Matrix {
   /* Compressed Sparse Column */
public:
   //CSCMatrix(size_t dimension, short fortran_indexing);
   CSCMatrix(size_t dimension, size_t maximum_number_nonzeros, short fortran_indexing);
   CSCMatrix(const std::vector<double>& matrix, const std::vector<size_t>& column_start, const std::vector<size_t>& row_number, size_t capacity,
         short fortran_indexing);

   std::vector<double> matrix;
   std::vector<size_t> column_start;
   std::vector<size_t> row_number;

   void insert(double term, size_t row_index, size_t column_index) override;
   std::vector<double> product(const std::vector<double>& vector) override;
   double quadratic_product(const std::vector<double>& x, const std::vector<double>& y) override;

   CSCMatrix add_identity_multiple(double multiple);
   double smallest_diagonal_entry();
   COOMatrix to_COO();
   UnoMatrix to_UnoMatrix(int uno_matrix_size);

   static CSCMatrix identity(size_t dimension, short fortran_indexing);

   friend std::ostream& operator<<(std::ostream& stream, CSCMatrix& matrix);
   friend std::ostream& operator<<(std::ostream& stream, const CSCMatrix& matrix);
};

class UnoMatrix : public Matrix {
   /* Coordinate list */
public:
   UnoMatrix(size_t dimension, size_t number_nonzeros, short fortran_indexing);

   SparseVector matrix;

   void insert(double term, size_t row_index, size_t column_index) override;
   std::vector<double> product(const std::vector<double>& vector) override;
   void add_matrix(UnoMatrix& other_matrix, double factor);

   double norm_1();
   COOMatrix to_COO();
   COOMatrix to_COO(const std::unordered_map<size_t, size_t>& mask);
   //CSCMatrix to_CSC();
   //CSCMatrix to_CSC(const std::unordered_map<int, int>& mask);

   friend std::ostream& operator<<(std::ostream& stream, UnoMatrix& matrix);
};

#endif // MATRIX_H
