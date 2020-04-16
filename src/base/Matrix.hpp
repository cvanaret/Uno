#ifndef MATRIX_H
#define MATRIX_H

#include <ostream>
#include <map>
#include <vector>

class Matrix {
public:
    Matrix(int size);
    virtual ~Matrix();

    int size;

    virtual int number_nonzeros() = 0;
    virtual std::vector<double> product(std::vector<double>& vector) = 0;
    virtual double quadratic_product(std::vector<double>& x, std::vector<double>& y) = 0;
};

class COOMatrix : public Matrix {
    /* Coordinate list */
public:
    COOMatrix(int size);

    std::vector<double> matrix;
    std::vector<int> row_indices;
    std::vector<int> column_indices;

    int number_nonzeros();
    /* build the matrix incrementally */
    void add_term(double term, int row_index, int column_index);
    double norm_1();
    std::vector<double> product(std::vector<double>& vector);
    double quadratic_product(std::vector<double>& x, std::vector<double>& y);

    /*COOMatrix add_identity_multiple(double multiple);*/
		
    friend std::ostream& operator<<(std::ostream &stream, COOMatrix& matrix);
    friend std::ostream& operator<<(std::ostream &stream, const COOMatrix& matrix);
};

// forward declaration
class CSCMatrix;

class ArgonotMatrix : public Matrix {
    /* Coordinate list */
public:
    ArgonotMatrix(int size);

    std::map<int, double> matrix;

    int number_nonzeros();
    /* build the matrix incrementally */
    void add_term(double term, int row_index, int column_index);
    double norm_1();
    std::vector<double> product(std::vector<double>& vector);
    double quadratic_product(std::vector<double>& x, std::vector<double>& y);
    COOMatrix to_COO();
    CSCMatrix to_CSC();
    
    friend std::ostream& operator<<(std::ostream &stream, ArgonotMatrix& matrix);
};

class CSCMatrix : public Matrix {
    /* Compressed Sparse Column */
public:
    CSCMatrix();
    CSCMatrix(std::vector<double>& matrix, std::vector<int>& column_start, std::vector<int>& row_number);

    std::vector<double> matrix;
    std::vector<int> column_start;
    std::vector<int> row_number;

    int number_nonzeros();
    std::vector<double> product(std::vector<double>& vector);
    double quadratic_product(std::vector<double>& x, std::vector<double>& y);
    CSCMatrix add_identity_multiple(double multiple);
    COOMatrix to_COO();
    ArgonotMatrix to_ArgonotMatrix(int argonot_matrix_size);

    friend std::ostream& operator<<(std::ostream &stream, CSCMatrix& matrix);
    friend std::ostream& operator<<(std::ostream &stream, const CSCMatrix& matrix);
};

#endif // MATRIX_H
