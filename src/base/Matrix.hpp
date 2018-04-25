#ifndef MATRIX_H
#define MATRIX_H

#include <ostream>
#include <vector>

class Matrix {
	public:
		Matrix();
		Matrix(std::vector<double>& matrix, std::vector<int>& column_start, std::vector<int>& row_number);
		
		std::vector<double> product(std::vector<double>& vector);
		double quadratic_product(std::vector<double>& x, std::vector<double>& y);
		Matrix add_identity_multiple(double multiple);
	
		std::vector<double> matrix;
		std::vector<int> column_start;
		std::vector<int> row_number;
		unsigned int number_nonzeros;
		
		friend std::ostream& operator<< (std::ostream &stream, Matrix& matrix);
		friend std::ostream& operator<< (std::ostream &stream, const Matrix& matrix);
};

#endif // MATRIX_H
