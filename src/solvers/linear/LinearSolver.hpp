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
    virtual std::vector<double> solve(const COOMatrix& matrix, std::vector<double>& rhs) = 0;

    [[nodiscard]] virtual std::tuple<int, int, int> get_inertia() const = 0;
    [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
    [[nodiscard]] virtual bool matrix_is_singular() const = 0;
    [[nodiscard]] virtual int rank() const = 0;
};

#endif // LINEARSOLVER_H
