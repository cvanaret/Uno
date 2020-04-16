#include <exception>
#include "Matrix.hpp"
#include "Utils.hpp"

Matrix::Matrix(int size) : size(size) {
}

Matrix::~Matrix() {
}

/* 
 * Coordinate list
 * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
 */

COOMatrix::COOMatrix(int size) : Matrix(size) {
}

int COOMatrix::number_nonzeros() {
    return this->matrix.size();
}

void COOMatrix::add_term(double term, int row_index, int column_index) {
    /* TODO: check matrix size */
    this->matrix.push_back(term);
    this->row_indices.push_back(row_index);
    this->column_indices.push_back(column_index);
    return;
}

double COOMatrix::norm_1() {
    // compute maximum column index
    int number_columns = 0;
    for (unsigned int k = 0; k < this->column_indices.size(); k++) {
        number_columns = std::max(number_columns, 1 + this->column_indices[k]);
    }
    // read the matrix and fill in the column_vectors norm vector
    std::vector<double> column_vectors(number_columns);
    for (unsigned int k = 0; k < this->matrix.size(); k++) {
        int j = this->column_indices[k];
        column_vectors[j] += std::abs(this->matrix[k]);
    }
    // compute the maximal component of the column_vectors vector
    double norm = 0.;
    for (unsigned int j = 0; j < column_vectors.size(); j++) {
        norm = std::max(norm, column_vectors[j]);
    }
    return norm;
}

std::vector<double> COOMatrix::product(std::vector<double>& vector) {
    throw std::runtime_error("COOMatrix::product not implemented");
}

double COOMatrix::quadratic_product(std::vector<double>& x, std::vector<double>& y) {
    throw std::runtime_error("COOMatrix::quadratic_product not implemented");
}

std::ostream& operator<<(std::ostream &stream, COOMatrix& matrix) {
    for (unsigned int k = 0; k < matrix.matrix.size(); k++) {
        stream << "m(" << matrix.row_indices[k] << ", " << matrix.column_indices[k] << ") = " << matrix.matrix[k] << "\n";
    }
    return stream;
}

std::ostream& operator<<(std::ostream &stream, const COOMatrix& matrix) {
    for (unsigned int k = 0; k < matrix.matrix.size(); k++) {
        stream << "m(" << matrix.row_indices[k] << ", " << matrix.column_indices[k] << ") = " << matrix.matrix[k] << "\n";
    }
    return stream;
}

/* 
 * Compressed Sparse Column
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 */

CSCMatrix::CSCMatrix() : Matrix(0) {
}

CSCMatrix::CSCMatrix(std::vector<double>& matrix, std::vector<int>& column_start, std::vector<int>& row_number) :
Matrix(column_start.size() - 1), matrix(matrix), column_start(column_start), row_number(row_number) {
}

int CSCMatrix::number_nonzeros() {
    return this->matrix.size();
}

CSCMatrix CSCMatrix::add_identity_multiple(double multiple) {
    /* initialize the damped matrix */
    std::vector<double> damped_matrix;
    std::vector<int> damped_column_start;
    std::vector<int> damped_row_number;

    int current_number_nonzeros = 0;
    damped_column_start.push_back(current_number_nonzeros);

    /* go through the columns */
    for (int j = 0; j < this->size; j++) {
        bool diagonal_term_updated = false;

        for (int k = this->column_start[j]; k < this->column_start[j + 1]; k++) {
            /* compute row number */
            int i = this->row_number[k];

            if (i == j) { /* update diagonal term */
                damped_matrix.push_back(this->matrix[k] + multiple);
                diagonal_term_updated = true;
            }
            else { /* keep off-diagonal term */
                damped_matrix.push_back(this->matrix[k]);
            }
            damped_row_number.push_back(i);
            current_number_nonzeros++;
        }
        /* add diagonal term at end of column if not present */
        if (!diagonal_term_updated) {
            damped_matrix.push_back(multiple);
            damped_row_number.push_back(j);
            current_number_nonzeros++;
        }
        damped_column_start.push_back(current_number_nonzeros);
    }
    return CSCMatrix(damped_matrix, damped_column_start, damped_row_number);
}

/* product of symmetric (n, n) matrix with (n, 1) vector */
std::vector<double> CSCMatrix::product(std::vector<double>& vector) {
    int n = this->column_start.size() - 1;
    /* create (n, 1) result */
    std::vector<double> result(n); // = {0.}
    for (int i = 0; i < n; i++) {
        result[i] = 0.;
    }

    for (int j = 0; j < n; j++) {
        for (int k = this->column_start[j]; k < this->column_start[j + 1]; k++) {
            int i = this->row_number[k];
            result[i] += vector[j] * this->matrix[k];

            /* non diagonal terms of the lower triangle */
            if (i != j) {
                result[j] += vector[i] * this->matrix[k];
            }
        }
    }
    return result;
}

/* compute x^T H y */
double CSCMatrix::quadratic_product(std::vector<double>& x, std::vector<double>& y) {
    if (x.size() != y.size()) {
        throw std::length_error("Matrix::quadratic_product: x and y have different sizes");
    }

    std::vector<double> hy = this->product(y); // H*y
    double product = dot(x, hy); // x^T*(H*y)
    return product;
}

COOMatrix CSCMatrix::to_COO() {
    COOMatrix coo_matrix(this->size);

    for (int j = 0; j < this->size; j++) {
        for (int k = this->column_start[j]; k < this->column_start[j + 1]; k++) {
            int i = this->row_number[k];
            coo_matrix.add_term(this->matrix[k], i, j);
        }
    }
    return coo_matrix;
}

ArgonotMatrix CSCMatrix::to_ArgonotMatrix(int argonot_matrix_size) {
    ArgonotMatrix argonot_matrix(argonot_matrix_size);
    
    for (int j = 0; j < this->size; j++) {
        for (int k = this->column_start[j]; k < this->column_start[j + 1]; k++) {
            int i = this->row_number[k];
            argonot_matrix.add_term(this->matrix[k], i, j);
        }
    }
    return argonot_matrix;
}

std::ostream& operator<<(std::ostream &stream, CSCMatrix& matrix) {
    /* Hessian */
    stream << "W = ";
    print_vector(stream, matrix.matrix, 0, 20);
    stream << "with column start: ";
    // TODO handle the stream
    print_vector(stream, matrix.column_start, 0, 20);
    stream << "and row number: ";
    print_vector(stream, matrix.row_number, 0, 20);
    return stream;
}

std::ostream& operator<<(std::ostream &stream, const CSCMatrix& matrix) {
    /* Hessian */
    stream << "W = ";
    print_vector(stream, matrix.matrix, 0, 20);
    stream << "with column start: ";
    // TODO handle the stream
    print_vector(stream, matrix.column_start, 0, 20);
    stream << "and row number: ";
    print_vector(stream, matrix.row_number, 0, 20);
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
 * Argonot matrix: bijection between indices and single key + sparse vector (map)
 */

ArgonotMatrix::ArgonotMatrix(int size) : Matrix(size) {
}

int ArgonotMatrix::number_nonzeros() {
    return this->matrix.size();
}

void ArgonotMatrix::add_term(double term, int row_index, int column_index) {
    // generate the unique key
    int key = column_index*this->size + row_index;
    // insert the element
    this->matrix[key] += term;
    return;
}

double ArgonotMatrix::norm_1() {
    // compute maximum column index
    int number_columns = 0;
    for (std::pair<const unsigned int, double> element: this->matrix) {
        int key = element.first;
        // retrieve indices
        int j = key/this->size;
        number_columns = std::max(number_columns, 1 + j);
    }
    // read the matrix and fill in the column_vectors norm vector
    std::vector<double> column_vectors(number_columns);
    for (std::pair<const unsigned int, double> element: this->matrix) {
        int key = element.first;
        double value = element.second;
        // retrieve indices
        int j = key/this->size;
        column_vectors[j] += std::abs(value);
    }
    // compute the maximal component of the column_vectors vector
    double norm = 0.;
    for (unsigned int j = 0; j < column_vectors.size(); j++) {
        norm = std::max(norm, column_vectors[j]);
    }
    return norm;
}

std::vector<double> ArgonotMatrix::product(std::vector<double>& vector) {
    std::vector<double> product;
    // TODO
    return product;
}

double ArgonotMatrix::quadratic_product(std::vector<double>& x, std::vector<double>& y) {
    // TODO
    return 0.;
}

COOMatrix ArgonotMatrix::to_COO() {
    COOMatrix coo_matrix(this->size);

    for (std::pair<const int, double> element: this->matrix) {
        int key = element.first;
        double value = element.second;
        // retrieve indices
        int i = key % this->size;
        int j = key / this->size;
        coo_matrix.add_term(value, i, j);
    }    
    return coo_matrix;
}

CSCMatrix ArgonotMatrix::to_CSC() {
    CSCMatrix csc_matrix;

    int previous_column = -1;
    int current_term = 0;
    for (std::pair<const int, double> element: this->matrix) {
        int key = element.first;
        double value = element.second;
        // retrieve indices
        int i = key % this->size;
        int j = key / this->size;
        
        csc_matrix.matrix.push_back(value);
        if (previous_column < j) {
            csc_matrix.column_start.push_back(current_term);
            previous_column = j;
        }
        csc_matrix.row_number.push_back(i);
        
        current_term++;
    }
    csc_matrix.column_start.push_back(current_term);
    return csc_matrix;
}

std::ostream& operator<<(std::ostream &stream, ArgonotMatrix& matrix) {
    for (std::pair<const int, double> element: matrix.matrix) {
        int key = element.first;
        double value = element.second;
        // retrieve indices
        int i = key % matrix.size;
        int j = key / matrix.size;
        stream << "m[" << i << ", " << j << "] = " << value << ", ";
    }
    return stream;
}

/*
void test_matrix() {
    int size = 3;
    int nnz = 4;
    ArgonotMatrix matrix(size, nnz);
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