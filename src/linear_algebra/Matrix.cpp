#include <exception>
#include <cmath>
#include <cassert>
#include <iomanip>
#include "Matrix.hpp"
#include "Vector.hpp"

Matrix::Matrix(size_t dimension, size_t capacity) : dimension(dimension), number_nonzeros(0), capacity(capacity) {
}

double Matrix::quadratic_product(const std::vector<double>& x, const std::vector<double>& y) {
   assert(x.size() == y.size() && "Matrix::quadratic_product: x and y have different sizes");

   std::vector<double> hy = this->product(y); // H*y
   double product = dot(x, hy); // x^T*(H*y)
   return product;
}

void Matrix::add_outer_product(const SparseVector& x, double scaling_factor) {
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
}

/* 
 * Coordinate list
 * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
 */

COOMatrix::COOMatrix(size_t dimension, size_t capacity) : Matrix(dimension, capacity) {
   this->matrix.reserve(capacity);
   this->row_indices.reserve(capacity);
   this->column_indices.reserve(capacity);
}

void COOMatrix::insert(double term, size_t row_index, size_t column_index) {
   /* TODO: check matrix size */
   this->matrix.push_back(term);
   this->row_indices.push_back((int) row_index);
   this->column_indices.push_back((int) column_index);
   this->number_nonzeros++;
}

double COOMatrix::norm_1() {
   // compute maximum column index
   size_t number_columns = 0;
   for (unsigned long j : this->column_indices) {
      number_columns = std::max(number_columns, 1 + j);
   }
   // read the matrix and fill in the column_vectors norm vector
   std::vector<double> column_vectors(number_columns);
   for (size_t k = 0; k < this->matrix.size(); k++) {
      size_t j = this->column_indices[k];
      column_vectors[j] += std::abs(this->matrix[k]);
   }
   // compute the maximal component of the column_vectors vector
   double norm = 0.;
   for (double & j : column_vectors) {
      norm = std::max(norm, j);
   }
   return norm;
}

// generic iterator
void COOMatrix::iter(const COOMatrix& matrix, const std::function<void (size_t, size_t, double)>& f) {
   for (size_t k = 0; k < matrix.number_nonzeros; k++) {
      size_t i = matrix.row_indices[k];
      size_t j = matrix.column_indices[k];
      f(i, j, matrix.matrix[k]);
   }
}

std::vector<double> COOMatrix::product(const std::vector<double>& vector) {
   std::vector<double> result(vector.size());
   COOMatrix::iter(*this, [&](size_t i, size_t j, double entry) {
      result[i] += entry * vector[j];
      // off-diagonal term
      if (i != j) {
         result[j] += entry * vector[i];
      }
   });
   return result;
}

std::ostream& operator<<(std::ostream& stream, COOMatrix& matrix) {
   for (size_t k = 0; k < matrix.number_nonzeros; k++) {
      stream << "m(" << matrix.row_indices[k] << ", " << matrix.column_indices[k] << ") = " << matrix.matrix[k] << "\n";
   }
   return stream;
}

std::ostream& operator<<(std::ostream& stream, const COOMatrix& matrix) {
   for (size_t k = 0; k < matrix.number_nonzeros; k++) {
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

CSCMatrix::CSCMatrix(size_t dimension, size_t capacity) : Matrix(dimension, capacity),
matrix(capacity), column_start(dimension + 1), row_number(capacity) {
}

CSCMatrix::CSCMatrix(const std::vector<double>& matrix, const std::vector<size_t>& column_start, const std::vector<size_t>& row_number, size_t
capacity) : Matrix(column_start.size() - 1, capacity), matrix(matrix), column_start(column_start), row_number(row_number) {
   assert(false && "CSCMatrix::CSCMatrix to check");
}

void CSCMatrix::insert(double /*term*/, size_t /*row_index*/, size_t /*column_index*/) {
   assert(false && "CSCMatrix::insert is not implemented");
}

// generic iterator
void CSCMatrix::iter(const CSCMatrix& matrix, const std::function<void (size_t, size_t, double)>& f) {
   for (size_t j = 0; j < matrix.dimension; j++) {
      for (size_t k = matrix.column_start[j]; k < matrix.column_start[j + 1]; k++) {
         size_t i = matrix.row_number[k];
         f(i, j, matrix.matrix[k]);
      }
   }
}

/* product of symmetric (n, n) matrix with (n, 1) vector */
std::vector<double> CSCMatrix::product(const std::vector<double>& vector) {
   assert(this->dimension == vector.size() && "the matrix and the vector do not have the same size");

   /* create (n, 1) result */
   std::vector<double> result(this->dimension); // = {0.}
   CSCMatrix::iter(*this, [&](size_t i, size_t j, double entry) {
      result[i] += vector[j] * entry;
      /* non diagonal terms of the lower triangle */
      if (i != j) {
         result[j] += vector[i] * entry;
      }
   });
   return result;
}

double CSCMatrix::quadratic_product(const std::vector<double>& x, const std::vector<double>& y) {
   assert(x.size() == y.size() && "CSCMatrix::quadratic_product: x and y have different sizes");

   double result = 0.;
   CSCMatrix::iter(*this, [&](size_t i, size_t j, double entry) {
      result += (i == j ? 1 : 2)*entry*x[i]*y[j];
   });
   return result;
}

CSCMatrix CSCMatrix::add_identity_multiple(double multiple) {
   /* initialize the damped matrix */
   std::vector<double> damped_matrix;
   std::vector<size_t> damped_column_start;
   std::vector<size_t> damped_row_number;
   damped_matrix.reserve(this->capacity);
   damped_row_number.reserve(this->capacity);
   damped_column_start.reserve(this->dimension);

   int current_number_nonzeros = 0;
   damped_column_start.push_back(0);

   /* go through the columns */
   for (size_t j = 0; j < this->dimension; j++) {
      bool diagonal_term_updated = false;

      for (size_t k = this->column_start[j]; k < this->column_start[j + 1]; k++) {
         /* compute row number */
         size_t i = this->row_number[k];

         if (i == j) { /* update diagonal term */
            damped_matrix.push_back(this->matrix[k] + multiple);
            diagonal_term_updated = true;
         }
         else if (j < i && !diagonal_term_updated) { /* we passed the diagonal (j, j) */
            damped_matrix.push_back(multiple);
            damped_matrix.push_back(this->matrix[k]);
            damped_row_number.push_back(j); // diagonal term
            current_number_nonzeros++;
            diagonal_term_updated = true;
         }
         else { /* keep off-diagonal term */
            damped_matrix.push_back(this->matrix[k]);
         }
         damped_row_number.push_back(i);
         current_number_nonzeros++;
      }
      /* add diagonal term in column j if not present */
      if (!diagonal_term_updated) {
         damped_matrix.push_back(multiple);
         damped_row_number.push_back(j);
         current_number_nonzeros++;
      }
      damped_column_start.push_back(current_number_nonzeros);
   }
   return CSCMatrix(damped_matrix, damped_column_start, damped_row_number, this->capacity);
}

double CSCMatrix::smallest_diagonal_entry() const {
   double smallest_entry = INFINITY;
   CSCMatrix::iter(*this, [&](size_t i, size_t j, double entry) {
      if (i == j) {
         smallest_entry = std::min(smallest_entry, entry);
      }
   });
   if (smallest_entry == INFINITY) {
      smallest_entry = 0.;
   }
   return smallest_entry;
}

COOMatrix CSCMatrix::to_COO() {
   COOMatrix coo_matrix(this->dimension, this->capacity);
   CSCMatrix::iter(*this, [&](size_t i, size_t j, double entry) {
      coo_matrix.insert(entry, i, j);
   });
   return coo_matrix;
}

UnoMatrix CSCMatrix::to_UnoMatrix(int uno_matrix_size) {
   UnoMatrix uno_matrix(uno_matrix_size, this->capacity);
   CSCMatrix::iter(*this, [&](size_t i, size_t j, double entry) {
      uno_matrix.insert(entry, i, j);
   });
   return uno_matrix;
}

CSCMatrix CSCMatrix::identity(size_t dimension) {
   /* initialize the identity matrix */
   std::vector<double> matrix(dimension);
   std::vector<size_t> column_start(dimension + 1);
   std::vector<size_t> row_number(dimension);

   for (size_t i = 0; i < dimension; i++) {
      matrix[i] = 1.;
      row_number[i] = i;
      column_start[i] = i;
   }
   column_start[dimension] = dimension;

   return CSCMatrix(matrix, column_start, row_number, dimension);
}

std::ostream& operator<<(std::ostream& stream, CSCMatrix& matrix) {
   stream << matrix.number_nonzeros << " non zeros\n";
   stream << "W = ";
   print_vector(stream, matrix.matrix, '\n', 0, matrix.number_nonzeros);
   stream << "with column start: ";
   print_vector(stream, matrix.column_start, '\n', 0, matrix.dimension + 1);
   stream << "and row number: ";
   print_vector(stream, matrix.row_number, '\n', 0, matrix.number_nonzeros);
   return stream;
}

std::ostream& operator<<(std::ostream& stream, const CSCMatrix& matrix) {
   stream << matrix.number_nonzeros << " non zeros\n";
   stream << "W = ";
   print_vector(stream, matrix.matrix, '\n', 0, matrix.number_nonzeros);
   stream << "with column start: ";
   print_vector(stream, matrix.column_start, '\n', 0, matrix.dimension + 1);
   stream << "and row number: ";
   print_vector(stream, matrix.row_number, '\n', 0, matrix.number_nonzeros);
   return stream;
}

/* 
 * Uno matrix: bijection between indices and single key + sparse vector (map)
 */

UnoMatrix::UnoMatrix(size_t dimension, size_t capacity): Matrix(dimension, capacity) {
}

void UnoMatrix::insert(double term, size_t row_index, size_t column_index) {
   // generate the unique key
   size_t key = column_index * this->dimension + row_index;
   // insert the element
   this->matrix[key] += term;
   this->number_nonzeros++;
}

double UnoMatrix::norm_1() {
   // compute maximum column index
   size_t number_columns = 0;
   for (std::pair<const size_t, double> element: this->matrix) {
      size_t key = element.first;
      // retrieve indices
      size_t j = key / this->dimension;
      number_columns = std::max(number_columns, 1 + j);
   }
   // read the matrix and fill in the column_vectors norm vector
   std::vector<double> column_vectors(number_columns);
   for (const auto[key, value]: this->matrix) {
      // retrieve indices
      size_t j = key / this->dimension;
      column_vectors[j] += std::abs(value);
   }
   // compute the maximal component of the column_vectors vector
   double norm = 0.;
   for (double column_vector : column_vectors) {
      norm = std::max(norm, column_vector);
   }
   return norm;
}

std::vector<double> UnoMatrix::product(const std::vector<double>& /*vector*/) {
   throw std::out_of_range("UnoMatrix::product is not implemented");
}

void UnoMatrix::add_matrix(UnoMatrix& other_matrix, double factor) {
   for (const auto[key, value]: other_matrix.matrix) {
      // retrieve indices
      size_t i = key % this->dimension;
      size_t j = key / this->dimension;
      this->insert(factor * value, i, j);
   }
}

COOMatrix UnoMatrix::to_COO() {
   COOMatrix coo_matrix(this->dimension, this->capacity);

   for (const auto[key, value]: this->matrix) {
      // retrieve indices
      size_t i = key % this->dimension;
      size_t j = key / this->dimension;
      coo_matrix.insert(value, i, j);
   }
   return coo_matrix;
}

/* generate a COO matrix by removing some variables (e.g. reduced Hessian in EQP problems) */
/* mask contains (i_origin, i_reduced) pairs, where i_origin is the original index, and i_reduced is the index in the reduced matrix */
COOMatrix UnoMatrix::to_COO(const std::unordered_map<size_t, size_t>& mask) {
   COOMatrix coo_matrix(this->dimension, this->capacity);

   for (const auto[key, value]: this->matrix) {
      // retrieve indices
      size_t i = key % this->dimension;
      size_t j = key / this->dimension;
      try {
         // if i and j are kept, compute the indices in the reduced matrix
         size_t i_reduced = mask.at(i);
         size_t j_reduced = mask.at(j);
         coo_matrix.insert(value, i_reduced, j_reduced);
      }
      catch (const std::out_of_range& e) {
      }
   }
   return coo_matrix;
}

//CSCMatrix UnoMatrix::to_CSC() {
//   CSCMatrix csc_matrix(this->dimension);
//
//   int current_column = 0;
//   int number_terms = 0;
//   csc_matrix.column_start.push_back(0);
//   for (const auto[key, value]: this->matrix) {
//      csc_matrix.matrix.push_back(value);
//      // retrieve indices
//      int i = key % this->dimension;
//      int j = key / this->dimension;
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
//   CSCMatrix csc_matrix(this->dimension);
//
//   int current_column = 0;
//   int number_terms = 0;
//   csc_matrix.column_start.push_back(0);
//   for (const auto[key, value]: this->matrix) {
//      // retrieve indices
//      int i = key % this->dimension;
//      int j = key / this->dimension;
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
      size_t i = key % matrix.dimension;
      size_t j = key / matrix.dimension;
      stream << "m(" << i << ", " << j << ") = " << value << ", ";
   }
   return stream;
}