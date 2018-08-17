#ifndef MATRIX_H
#define MATRIX_H

#include <ostream>
#include <vector>

class Matrix {
	public:
		Matrix(unsigned int size, unsigned int number_nonzeros);
		virtual ~Matrix();
		
		unsigned int size;
		unsigned int number_nonzeros;
		
		virtual std::vector<double> product(std::vector<double>& vector) = 0;
		virtual double quadratic_product(std::vector<double>& x, std::vector<double>& y) = 0;
};

class COOMatrix: public Matrix {
	/* Coordinate list */
	public:
		COOMatrix(unsigned int size, unsigned int number_nonzeros);
		
		std::vector<double> matrix;
		std::vector<int> row_indices;
		std::vector<int> column_indices;
		
		/* build the matrix incrementally */
		void add_term(double term, int row_index, int column_index);
		std::vector<double> product(std::vector<double>& vector);
		double quadratic_product(std::vector<double>& x, std::vector<double>& y);
		
		/*COOMatrix add_identity_multiple(double multiple);
		
		friend std::ostream& operator<< (std::ostream &stream, COOMatrix& matrix);
		friend std::ostream& operator<< (std::ostream &stream, const COOMatrix& matrix);*/
};

class CSCMatrix: public Matrix {
	/* Compressed Sparse Column */
	public:
		CSCMatrix();
		CSCMatrix(std::vector<double>& matrix, std::vector<int>& column_start, std::vector<int>& row_number);
		
		std::vector<double> matrix;
		std::vector<int> column_start;
		std::vector<int> row_number;
		
		std::vector<double> product(std::vector<double>& vector);
		double quadratic_product(std::vector<double>& x, std::vector<double>& y);
		CSCMatrix add_identity_multiple(double multiple);
		COOMatrix to_COO();
		
		friend std::ostream& operator<< (std::ostream &stream, CSCMatrix& matrix);
		friend std::ostream& operator<< (std::ostream &stream, const CSCMatrix& matrix);
};

#endif // MATRIX_H
