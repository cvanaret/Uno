#ifndef MATRIX_H
#define MATRIX_H

#include <ostream>
#include <map>
#include <vector>

class Matrix {
public:
    Matrix(int dimension, int fortran_indexing);
    virtual ~Matrix();

    int dimension;
    int fortran_indexing;

    virtual int number_nonzeros() = 0;
    /* build the matrix incrementally */
    virtual void add_term(double term, int row_index, int column_index) = 0;
    virtual std::vector<double> product(std::vector<double>& vector) = 0;
    
    double quadratic_product(std::vector<double>& x, std::vector<double>& y);
    void add_outer_product(std::map<int, double>& x, double scaling_factor = 1.);
};

class COOMatrix : public Matrix {
    /* Coordinate list */
public:
    COOMatrix(int dimension, int fortran_indexing);

    std::vector<double> matrix;
    std::vector<int> row_indices;
    std::vector<int> column_indices;

    int number_nonzeros() override;
    void add_term(double term, int row_index, int column_index) override;
    std::vector<double> product(std::vector<double>& vector) override;
    
    double norm_1();
    /*COOMatrix add_identity_multiple(double multiple);*/

    friend std::ostream& operator<<(std::ostream &stream, COOMatrix& matrix);
    friend std::ostream& operator<<(std::ostream &stream, const COOMatrix& matrix);
};

// forward declaration
class ArgonotMatrix;

class CSCMatrix : public Matrix {
    /* Compressed Sparse Column */
public:
    CSCMatrix();
    CSCMatrix(std::vector<double>& matrix, std::vector<int>& column_start, std::vector<int>& row_number, int fortran_indexing);

    std::vector<double> matrix;
    std::vector<int> column_start;
    std::vector<int> row_number;

    int number_nonzeros() override;
    void add_term(double term, int row_index, int column_index) override;
    std::vector<double> product(std::vector<double>& vector) override;
    
    CSCMatrix add_identity_multiple(double multiple);
    double smallest_diagonal_entry();
    COOMatrix to_COO();
    ArgonotMatrix to_ArgonotMatrix(int argonot_matrix_size);
    
    static CSCMatrix identity(int dimension, int fortran_indexing);

    friend std::ostream& operator<<(std::ostream &stream, CSCMatrix& matrix);
    friend std::ostream& operator<<(std::ostream &stream, const CSCMatrix& matrix);
};

class ArgonotMatrix : public Matrix {
    /* Coordinate list */
public:
    ArgonotMatrix(int dimension, int fortran_indexing);

    std::map<int, double> matrix;

    int number_nonzeros() override;
    void add_term(double term, int row_index, int column_index) override;
    std::vector<double> product(std::vector<double>& vector) override;
    
    double norm_1();
    COOMatrix to_COO();
    CSCMatrix to_CSC();

    friend std::ostream& operator<<(std::ostream &stream, ArgonotMatrix& matrix);
};

#endif // MATRIX_H
