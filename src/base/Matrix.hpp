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

    virtual int number_nonzeros() = 0;
    /* build the matrix incrementally */
    virtual void insert(double term, int row_index, int column_index) = 0;
    virtual std::vector<double> product(std::vector<double>& vector) = 0;
    
    double quadratic_product(std::vector<double>& x, std::vector<double>& y);
    void add_outer_product(SparseGradient& x, double scaling_factor = 1.);
};

class COOMatrix : public Matrix {
    /* Coordinate list */
public:
    COOMatrix(int dimension, short fortran_indexing);

    std::vector<double> matrix;
    std::vector<int> row_indices;
    std::vector<int> column_indices;

    int number_nonzeros() override;
    void insert(double term, int row_index, int column_index) override;
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
    CSCMatrix(int dimension, short fortran_indexing);
    CSCMatrix(std::vector<double>& matrix, std::vector<int>& column_start, std::vector<int>& row_number, int fortran_indexing);

    std::vector<double> matrix;
    std::vector<int> column_start;
    std::vector<int> row_number;

    int number_nonzeros() override;
    void insert(double term, int row_index, int column_index) override;
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
    ArgonotMatrix(int dimension, short fortran_indexing);

    SparseGradient matrix;

    int number_nonzeros() override;
    void insert(double term, int row_index, int column_index) override;
    std::vector<double> product(std::vector<double>& vector) override;
    void add_matrix(ArgonotMatrix& other_matrix, double factor);
    
    double norm_1();
    COOMatrix to_COO();
    COOMatrix to_COO(std::unordered_map<int, int> mask);
    CSCMatrix to_CSC();
    CSCMatrix to_CSC(std::unordered_map<int, int> mask);

    friend std::ostream& operator<<(std::ostream &stream, ArgonotMatrix& matrix);
};

#endif // MATRIX_H
