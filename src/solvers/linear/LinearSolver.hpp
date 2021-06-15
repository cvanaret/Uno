#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <vector>
#include "Matrix.hpp"

class LinearSolver {
public:

    LinearSolver() = default;
    virtual ~LinearSolver() = default;
    virtual void factorize(COOMatrix& matrix) = 0;
    virtual void do_symbolic_factorization(const COOMatrix& matrix) = 0;
    virtual void do_numerical_factorization(const COOMatrix& matrix) = 0;
    virtual void solve(std::vector<double>& rhs) = 0;
    
    virtual size_t number_negative_eigenvalues() const = 0;
    virtual bool matrix_is_singular() const = 0;
    virtual int rank() const = 0;
};

#endif // LINEARSOLVER_H
