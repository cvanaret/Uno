#ifndef MATRIX_H
#define MATRIX_H

#include <ostream>
#include <map>
#include <vector>
#include "SparseGradient.hpp"

class Matrix {
public:
   Matrix(int dimension, short fortran_indexing);
   virtual ~Matrix();

   int dimension;
   short fortran_indexing;

   virtual size_t number_nonzeros() const = 0;
   /* build the matrix incrementally */
   virtual void insert(double term, int row_index, int column_index) = 0;
   virtual std::vector<double> product(const std::vector<double>& vector) = 0;

   double quadratic_product(const std::vector<double>& x, const std::vector<double>& y);
   void add_outer_product(const SparseGradient& x, double scaling_factor = 1.);
};

class COOMatrix : public Matrix {
   /* Coordinate list */
public:
   COOMatrix(int dimension, short fortran_indexing);

   std::vector<double> matrix;
   std::vector<int> row_indices;
   std::vector<int> column_indices;

   size_t number_nonzeros() const override;
   void insert(double term, int row_index, int column_index) override;
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
   CSCMatrix(int dimension, short fortran_indexing);
   CSCMatrix(const std::vector<double>& matrix, const std::vector<int>& column_start, const std::vector<int>& row_number,
         int fortran_indexing);

   std::vector<double> matrix;
   std::vector<int> column_start;
   std::vector<int> row_number;

   size_t number_nonzeros() const override;
   void insert(double term, int row_index, int column_index) override;
   std::vector<double> product(const std::vector<double>& vector) override;

   CSCMatrix add_identity_multiple(double multiple);
   double smallest_diagonal_entry();
   COOMatrix to_COO();
   UnoMatrix to_UnoMatrix(int uno_matrix_size);

   static CSCMatrix identity(int dimension, int fortran_indexing);

   friend std::ostream& operator<<(std::ostream& stream, CSCMatrix& matrix);
   friend std::ostream& operator<<(std::ostream& stream, const CSCMatrix& matrix);
};

class UnoMatrix : public Matrix {
   /* Coordinate list */
public:
   UnoMatrix(int dimension, short fortran_indexing);

   SparseGradient matrix;

   size_t number_nonzeros() const override;
   void insert(double term, int row_index, int column_index) override;
   std::vector<double> product(const std::vector<double>& vector) override;
   void add_matrix(UnoMatrix& other_matrix, double factor);

   double norm_1();
   COOMatrix to_COO();
   COOMatrix to_COO(const std::unordered_map<int, int>& mask);
   CSCMatrix to_CSC();
   CSCMatrix to_CSC(const std::unordered_map<int, int>& mask);

   friend std::ostream& operator<<(std::ostream& stream, UnoMatrix& matrix);
};

#endif // MATRIX_H
