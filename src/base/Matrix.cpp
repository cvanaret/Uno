#include "Matrix.hpp"
#include "Utils.hpp"

Matrix::Matrix(unsigned int size, unsigned int number_nonzeros): size(size), number_nonzeros(number_nonzeros) {
}

Matrix::~Matrix() {
}

/* 
 * Coordinate list
 * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
*/

COOMatrix::COOMatrix(unsigned int size, unsigned int number_nonzeros): Matrix(size, number_nonzeros) {
}

void COOMatrix::add_term(double term, int row_index, int column_index) {
	/* TODO: check matrix size */
	this->matrix.push_back(term);
	this->row_indices.push_back(row_index);
	this->column_indices.push_back(column_index);
	this->number_nonzeros++;
	
	return;
}

std::vector<double> COOMatrix::product(std::vector<double>& vector) {
	throw std::invalid_argument("COOMatrix::product not implemented");
}

double COOMatrix::quadratic_product(std::vector<double>& x, std::vector<double>& y) {
	throw std::invalid_argument("COOMatrix::quadratic_product not implemented");
}

/* 
 * Compressed Sparse Column
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
*/

CSCMatrix::CSCMatrix(): Matrix(0, 0) {
}

CSCMatrix::CSCMatrix(std::vector<double>& matrix, std::vector<int>& column_start, std::vector<int>& row_number):
		Matrix(column_start.size()-1, matrix.size()), matrix(matrix), column_start(column_start), row_number(row_number) {
}

CSCMatrix CSCMatrix::add_identity_multiple(double multiple) {
	/* initialize the damped matrix */
	std::vector<double> damped_matrix;
	std::vector<int> damped_column_start;
	std::vector<int> damped_row_number;
	
	int current_number_nonzeros = 0;
	damped_column_start.push_back(current_number_nonzeros);
	
	/* go through the columns */
	for (unsigned int j = 0; j < this->size; j++) {
		bool diagonal_term_updated = false;
		
		for (int k = this->column_start[j]; k < this->column_start[j+1]; k++) {
			/* compute row number */
			unsigned int i = this->row_number[k];
			
			if (i == j) { /* update diagonal term */
				damped_matrix.push_back(this->matrix[k] + multiple);
				diagonal_term_updated = true;
			}
			else { /* keep off-diagonal term */
				damped_matrix.push_back(this->matrix[k]);
			}
			damped_row_number.push_back(this->row_number[k]);
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
	int n = this->column_start.size()-1;
	/* create (n, 1) result */
	std::vector<double> result(n); // = {0.}
	for (int i = 0; i < n; i++) {
		result[i] = 0.;
	}
	
	for (int j = 0; j < n; j++) {
		for (int k = this->column_start[j]; k < this->column_start[j+1]; k++) {
			int i = this->row_number[k];
			result[i] += vector[j]*this->matrix[k];
			
			/* non diagonal terms of the lower triangle */
			if (i != j) {
				result[j] += vector[i]*this->matrix[k];
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
	COOMatrix coo_matrix(this->size, this->number_nonzeros);
	
	for (unsigned int j = 0; j < this->size; j++) {
		for (int k = this->column_start[j]; k < this->column_start[j+1]; k++) {
			int i = this->row_number[k];
			
			/* TODO */
			coo_matrix.row_indices.push_back(i);
			coo_matrix.column_indices.push_back(j);
			coo_matrix.matrix.push_back(this->matrix[k]);
		}
	}
	return coo_matrix;
}

std::ostream& operator<< (std::ostream &stream, CSCMatrix& matrix) {
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

std::ostream& operator<< (std::ostream &stream, const CSCMatrix& matrix) {
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
