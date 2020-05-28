#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <vector>
#include "Matrix.hpp"

class LinearSolver {
public:

    virtual ~LinearSolver() {
    };
    virtual void factorize(COOMatrix& matrix) = 0;
    virtual void solve(std::vector<double>& rhs) = 0;
    virtual int number_negative_eigenvalues() = 0;
    virtual bool matrix_is_singular() = 0;
    virtual int rank() = 0;
};

#endif // LINEARSOLVER_H
