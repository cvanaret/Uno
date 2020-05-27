#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <vector>
#include "Matrix.hpp"

struct MA57Factorization {
    int dimension;
	std::vector<double> fact;
	int lfact;
	std::vector<int> ifact;
	int lifact;
	std::vector<int> iwork;
    std::vector<int> info;
    
    int number_negative_eigenvalues();
    bool matrix_is_singular();
    int rank();
};

class LinearSolver {
public:

    virtual ~LinearSolver() {
    };
    virtual MA57Factorization factorize(COOMatrix& matrix) = 0;
    virtual void solve(MA57Factorization& factorization, std::vector<double>& rhs) = 0;
    virtual void solve(COOMatrix& matrix, std::vector<double>& rhs) = 0;
};

#endif // LINEARSOLVER_H
