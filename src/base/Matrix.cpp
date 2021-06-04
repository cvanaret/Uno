#include <exception>
#include <cmath>
#include "Matrix.hpp"
#include "Vector.hpp"

Matrix::Matrix(size_t dimension, size_t number_nonzeros, short fortran_indexing) :
dimension(dimension), number_nonzeros(number_nonzeros), fortran_indexing(fortran_indexing) {
}

Matrix::~Matrix() {
}

double Matrix::quadratic_product(const std::vector<double>& x, const std::vector<double>& y) {
   if (x.size() != y.size()) {
      throw std::length_error("COOMatrix::quadratic_product: x and y have different sizes");
   }

   std::vector<double> hy = this->product(y); // H*y
   double product = dot(x, hy); // x^T*(H*y)
   return product;
}

void Matrix::add_outer_product(const SparseGradient& x, double scaling_factor) {
   /* perform matrix addition: A + rho x x^T */
   for (const auto[row_index, row_term]: x) {
      for (const auto[column_index, column_term]: x) {
         // upper triangular matrix
         if (row_index <= column_index) {
            // add product of components
            this->insert(scaling_factor * row_term * column_term, row_index, column_index);
         }
      }
   }
   return;
}

/* 
 * Coordinate list
 * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
 */

COOMatrix::COOMatrix(size_t dimension, size_t number_nonzeros, short fortran_indexing) :
Matrix(dimension, number_nonzeros, fortran_indexing) {
}

void COOMatrix::insert(double term, int row_index, int column_index) {
   /* TODO: check matrix size */
   this->matrix.push_back(term);
   this->row_indices.push_back(row_index + this->fortran_indexing);
   this->column_indices.push_back(column_index + this->fortran_indexing);
   this->number_nonzeros++;
   return;
}

double COOMatrix::norm_1() {
   // compute maximum column index
   int number_columns = 0;
   for (size_t k = 0; k < this->column_indices.size(); k++) {
      number_columns = std::max(number_columns, 1 + this->column_indices[k] - this->fortran_indexing);
   }
   // read the matrix and fill in the column_vectors norm vector
   std::vector<double> column_vectors(number_columns);
   for (size_t k = 0; k < this->matrix.size(); k++) {
      int j = this->column_indices[k] - this->fortran_indexing;
      column_vectors[j] += std::abs(this->matrix[k]);
   }
   // compute the maximal component of the column_vectors vector
   double norm = 0.;
   for (size_t j = 0; j < column_vectors.size(); j++) {
      norm = std::max(norm, column_vectors[j]);
   }
   return norm;
}

std::vector<double> COOMatrix::product(const std::vector<double>& vector) {
   std::vector<double> result(vector.size());
   for (size_t k = 0; k < this->matrix.size(); k++) {
      int i = this->row_indices[k] - this->fortran_indexing;
      int j = this->column_indices[k] - this->fortran_indexing;
      result[i] += this->matrix[k] * vector[j];

      // off-diagonal term
      if (i != j) {
         result[j] += this->matrix[k] * vector[i];
      }
   }
   return result;
}

std::ostream& operator<<(std::ostream& stream, COOMatrix& matrix) {
   for (size_t k = 0; k < matrix.matrix.size(); k++) {
      stream << "m(" << matrix.row_indices[k] << ", " << matrix.column_indices[k] << ") = " << matrix.matrix[k] << "\n";
   }
   return stream;
}

std::ostream& operator<<(std::ostream& stream, const COOMatrix& matrix) {
   for (size_t k = 0; k < matrix.matrix.size(); k++) {
      stream << "m(" << matrix.row_indices[k] << ", " << matrix.column_indices[k] << ") = " << matrix.matrix[k] << "\n";
   }
   return stream;
}

/* 
 * Compressed Sparse Column
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 */

// matrix and row_number have nnz elements
// column_start has dimension+1 elements

CSCMatrix::CSCMatrix(size_t dimension, int maximum_number_nonzeros, int fortran_indexing) :
Matrix(dimension, maximum_number_nonzeros, fortran_indexing),
matrix(maximum_number_nonzeros), column_start(dimension + 1), row_number(maximum_number_nonzeros) {
}

CSCMatrix::CSCMatrix(const std::vector<double>& matrix, const std::vector<int>& column_start, const std::vector<int>& row_number,
      int fortran_indexing) :
      Matrix(column_start.size() - 1, matrix.size(), fortran_indexing),
      matrix(matrix), column_start(column_start), row_number(row_number) {
}

void CSCMatrix::insert(double /*term*/, int /*row_index*/, int /*column_index*/) {
   throw std::out_of_range("CSCMatrix::add_term is not implemented");
}

/* product of symmetric (n, n) matrix with (n, 1) vector */
std::vector<double> CSCMatrix::product(const std::vector<double>& vector) {
   size_t n = this->column_start.size() - 1;
   /* create (n, 1) result */
   std::vector<double> result(n); // = {0.}

   for (size_t j = 0; j < n; j++) {
      for (int k = this->column_start[j] - this->fortran_indexing; k < this->column_start[j + 1] - this->fortran_indexing; k++) {
         size_t i = this->row_number[k] - this->fortran_indexing;
         result[i] += vector[j] * this->matrix[k];

         /* non diagonal terms of the lower triangle */
         if (i != j) {
            result[j] += vector[i] * this->matrix[k];
         }
      }
   }
   return result;
}

CSCMatrix CSCMatrix::add_identity_multiple(double multiple) {
   /* initialize the damped matrix */
   std::vector<double> damped_matrix;
   std::vector<int> damped_column_start;
   std::vector<int> damped_row_number;

   int current_number_nonzeros = 0;
   damped_column_start.push_back(this->fortran_indexing);

   /* go through the columns */
   for (size_t j = 0; j < this->dimension; j++) {
      bool diagonal_term_updated = false;

      for (int k = this->column_start[j] - this->fortran_indexing; k < this->column_start[j + 1] - this->fortran_indexing; k++) {
         /* compute row number */
         size_t i = this->row_number[k] - this->fortran_indexing;

         if (i == j) { /* update diagonal term */
            damped_matrix.push_back(this->matrix[k] + multiple);
            diagonal_term_updated = true;
         }
         else if (j < i && !diagonal_term_updated) { /* we passed the diagonal (j, j) */
            damped_matrix.push_back(multiple);
            damped_matrix.push_back(this->matrix[k]);
            damped_row_number.push_back(j + this->fortran_indexing); // diagonal term
            current_number_nonzeros++;
            diagonal_term_updated = true;
         }
         else { /* keep off-diagonal term */
            damped_matrix.push_back(this->matrix[k]);
         }
         damped_row_number.push_back(i + this->fortran_indexing);
         current_number_nonzeros++;
      }
      /* add diagonal term in column j if not present */
      if (!diagonal_term_updated) {
         damped_matrix.push_back(multiple);
         damped_row_number.push_back(j + this->fortran_indexing);
         current_number_nonzeros++;
      }
      damped_column_start.push_back(current_number_nonzeros + this->fortran_indexing);
   }
   return CSCMatrix(damped_matrix, damped_column_start, damped_row_number, this->fortran_indexing);
}

double CSCMatrix::smallest_diagonal_entry() {
   double smallest_entry = INFINITY;

   for (size_t j = 0; j < this->dimension; j++) {
      for (int k = this->column_start[j] - this->fortran_indexing; k < this->column_start[j + 1] - this->fortran_indexing; k++) {
         /* compute row number */
         size_t i = this->row_number[k] - this->fortran_indexing;

         if (i == j) {
            smallest_entry = std::min(smallest_entry, this->matrix[k]);
         }
      }
   }

   if (smallest_entry == INFINITY) {
      smallest_entry = 0.;
   }
   return smallest_entry;
}

COOMatrix CSCMatrix::to_COO() {
   COOMatrix coo_matrix(this->dimension, this->number_nonzeros, this->fortran_indexing);

   for (size_t j = 0; j < this->dimension; j++) {
      for (int k = this->column_start[j] - this->fortran_indexing; k < this->column_start[j + 1] - this->fortran_indexing; k++) {
         int i = this->row_number[k] - this->fortran_indexing;
         coo_matrix.insert(this->matrix[k], i, j);
      }
   }
   return coo_matrix;
}

UnoMatrix CSCMatrix::to_UnoMatrix(int uno_matrix_size) {
   UnoMatrix uno_matrix(uno_matrix_size, this->number_nonzeros, this->fortran_indexing);

   for (size_t j = 0; j < this->dimension; j++) {
      for (int k = this->column_start[j] - this->fortran_indexing; k < this->column_start[j + 1] - this->fortran_indexing; k++) {
         int i = this->row_number[k] - this->fortran_indexing;
         uno_matrix.insert(this->matrix[k], i, j);
      }
   }
   return uno_matrix;
}

CSCMatrix CSCMatrix::identity(size_t dimension, int fortran_indexing) {
   /* initialize the identity matrix */
   std::vector<double> matrix(dimension);
   std::vector<int> column_start(dimension + 1);
   std::vector<int> row_number(dimension);

   for (size_t i = 0; i < dimension; i++) {
      matrix[i] = 1.;
      row_number[i] = i + fortran_indexing;
      column_start[i] = i + fortran_indexing;
   }
   column_start[dimension] = dimension + fortran_indexing;

   return CSCMatrix(matrix, column_start, row_number, fortran_indexing);
}

std::ostream& operator<<(std::ostream& stream, CSCMatrix& matrix) {
   /* Hessian */
   stream << "W = ";
   print_vector(stream, matrix.matrix);
   stream << "with column start: ";
   print_vector(stream, matrix.column_start);
   stream << "and row number: ";
   print_vector(stream, matrix.row_number);
   return stream;
}

std::ostream& operator<<(std::ostream& stream, const CSCMatrix& matrix) {
   /* Hessian */
   stream << "W = ";
   print_vector(stream, matrix.matrix);
   stream << "with column start: ";
   print_vector(stream, matrix.column_start);
   stream << "and row number: ";
   print_vector(stream, matrix.row_number);
   return stream;
}

//// test of add_identity_multiple
//std::vector<double> matrix(6);
//matrix[0] = matrix[1] = matrix[2] = matrix[3] = matrix[4] = matrix[5] = 1.;
//std::vector<int> column_start(5);
//column_start[0] = 1;
//column_start[1] = 2;
//column_start[2] = 3;
//column_start[3] = 4;
//column_start[4] = 7;
//std::vector<int> row_number(6);
//row_number[0] = 1;
//row_number[1] = 1;
//row_number[2] = 2;
//row_number[3] = 1;
//row_number[4] = 3;
//row_number[5] = 4;
//Matrix m(matrix, column_start, row_number);
//m.display();

//Matrix m_damped = add_identity_multiple(m, 100.);
//m_damped.display();

//// expect 
//W = 101 1 100 1 100 1 1 101 
//with column start: 1 2 4 6 9 
//and row number: 1 1 2 2 3 1 3 4 

//void test_product() {
//std::vector<double> m(5);
//m[0] = 1.;
//m[1] = 3.;
//m[2] = 4.;
//m[3] = 7.;
//m[4] = 15.;
//std::vector<int> matrix_column_start = {1, 2, 4, 6};
//std::vector<int> matrix_row_number = {1, 1, 2, 1, 3};
//Matrix matrix(m, matrix_column_start, matrix_row_number);
//std::vector<double> vector(3);
//vector[0] = 1.;
//vector[1] = 2.;
//vector[2] = 3.;

//std::vector<double> result = matrix.product(vector);
///* *expected result = (28, 11, 52) */
//}

/* 
 * Uno matrix: bijection between indices and single key + sparse vector (map)
 */

UnoMatrix::UnoMatrix(size_t dimension, size_t number_nonzeros, short fortran_indexing) :
   Matrix(dimension, number_nonzeros, fortran_indexing) {
}

void UnoMatrix::insert(double term, int row_index, int column_index) {
   // generate the unique key
   int key = column_index * this->dimension + row_index;
   // insert the element
   this->matrix[key] += term;
   return;
}

double UnoMatrix::norm_1() {
   // compute maximum column index
   int number_columns = 0;
   for (std::pair<const size_t, double> element: this->matrix) {
      int key = element.first;
      // retrieve indices
      int j = key / this->dimension;
      number_columns = std::max(number_columns, 1 + j);
   }
   // read the matrix and fill in the column_vectors norm vector
   std::vector<double> column_vectors(number_columns);
   for (const auto[key, value]: this->matrix) {
      // retrieve indices
      int j = key / this->dimension;
      column_vectors[j] += std::abs(value);
   }
   // compute the maximal component of the column_vectors vector
   double norm = 0.;
   for (size_t j = 0; j < column_vectors.size(); j++) {
      norm = std::max(norm, column_vectors[j]);
   }
   return norm;
}

std::vector<double> UnoMatrix::product(const std::vector<double>& /*vector*/) {
   throw std::out_of_range("UnoMatrix::product is not implemented");
}

void UnoMatrix::add_matrix(UnoMatrix& other_matrix, double factor) {
   for (const auto[key, value]: other_matrix.matrix) {
      // retrieve indices
      int i = key % this->dimension;
      int j = key / this->dimension;
      this->insert(factor * value, i, j);
   }
   return;
}

COOMatrix UnoMatrix::to_COO() {
   COOMatrix coo_matrix(this->dimension, this->number_nonzeros, this->fortran_indexing);

   for (const auto[key, value]: this->matrix) {
      // retrieve indices
      int i = key % this->dimension;
      int j = key / this->dimension;
      coo_matrix.insert(value, i, j);
   }
   return coo_matrix;
}

/* generate a COO matrix by removing some variables (e.g. reduced Hessian in EQP problems) */
/* mask contains (i_origin, i_reduced) pairs, where i_origin is the original index, and i_reduced is the index in the reduced matrix */
COOMatrix UnoMatrix::to_COO(const std::unordered_map<int, int>& mask) {
   COOMatrix coo_matrix(this->dimension, this->number_nonzeros, this->fortran_indexing);

   for (const auto[key, value]: this->matrix) {
      // retrieve indices
      int i = key % this->dimension;
      int j = key / this->dimension;
      try {
         // if i and j are kept, compute the indices in the reduced matrix
         int i_reduced = mask.at(i);
         int j_reduced = mask.at(j);
         coo_matrix.insert(value, i_reduced, j_reduced);
      }
      catch (const std::out_of_range& e) {
      }
   }
   return coo_matrix;
}

//CSCMatrix UnoMatrix::to_CSC() {
//   CSCMatrix csc_matrix(this->dimension, this->fortran_indexing);
//
//   int current_column = this->fortran_indexing;
//   int number_terms = this->fortran_indexing;
//   csc_matrix.column_start.push_back(this->fortran_indexing);
//   for (const auto[key, value]: this->matrix) {
//      csc_matrix.matrix.push_back(value);
//      // retrieve indices
//      int i = key % this->dimension + this->fortran_indexing;
//      int j = key / this->dimension + this->fortran_indexing;
//
//      csc_matrix.row_number.push_back(i);
//      for (int column = current_column; column < j; column++) {
//         csc_matrix.column_start.push_back(number_terms);
//      }
//      current_column = j;
//      number_terms++;
//   }
//   csc_matrix.column_start.push_back(number_terms);
//   return csc_matrix;
//}

//CSCMatrix UnoMatrix::to_CSC(const std::unordered_map<int, int>& mask) {
//   CSCMatrix csc_matrix(this->dimension, this->fortran_indexing);
//
//   int current_column = this->fortran_indexing;
//   int number_terms = this->fortran_indexing;
//   csc_matrix.column_start.push_back(this->fortran_indexing);
//   for (const auto[key, value]: this->matrix) {
//      // retrieve indices
//      int i = key % this->dimension + this->fortran_indexing;
//      int j = key / this->dimension + this->fortran_indexing;
//      try {
//         // if i and j are kept, compute the indices in the reduced matrix
//         int i_reduced = mask.at(i);
//         int j_reduced = mask.at(j);
//         csc_matrix.matrix.push_back(value);
//         csc_matrix.row_number.push_back(i_reduced);
//         for (int column = current_column; column < j_reduced; column++) {
//            csc_matrix.column_start.push_back(number_terms);
//         }
//         current_column = j_reduced;
//         number_terms++;
//      }
//      catch (const std::out_of_range& e) {
//      }
//   }
//   csc_matrix.column_start.push_back(number_terms);
//   return csc_matrix;
//}

std::ostream& operator<<(std::ostream& stream, UnoMatrix& matrix) {
   for (const auto[key, value]: matrix.matrix) {
      // retrieve indices
      int i = key % matrix.dimension;
      int j = key / matrix.dimension;
      stream << "m(" << i << ", " << j << ") = " << value << ", ";
   }
   return stream;
}

/*
void test_matrix() {
    int size = 3;
    int nnz = 4;
    UnoMatrix matrix(size, nnz);
    // add terms
    matrix.add_term(12., 0, 0);
    
    matrix.add_term(14., 0, 1);
    matrix.add_term(7., 0, 1);
    
    matrix.add_term(13., 1, 2);
    
    matrix.add_term(19., 2, 2);
    
    print_vector(std::cout, matrix.matrix);
    
    COOMatrix coo_matrix = matrix.to_COO();
    std::cout << "COO matrix: \n";
    print_vector(std::cout, coo_matrix.matrix);
    print_vector(std::cout, coo_matrix.row_indices);
    print_vector(std::cout, coo_matrix.column_indices);
    
    return;
}
 */
