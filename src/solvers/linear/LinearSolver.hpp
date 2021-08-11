#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <vector>
#include "COOSymmetricMatrix.hpp"
#include "CSCSymmetricMatrix.hpp"

class LinearSolver {
public:
   LinearSolver() = default;
   virtual ~LinearSolver() = default;
   virtual void factorize(COOSymmetricMatrix& matrix) = 0;
   virtual void do_symbolic_factorization(COOSymmetricMatrix& matrix) = 0;
   virtual void do_numerical_factorization(COOSymmetricMatrix& matrix) = 0;
   virtual std::vector<double> solve(COOSymmetricMatrix& matrix, const std::vector<double>& rhs) = 0;

   virtual void factorize(CSCSymmetricMatrix& matrix) = 0;
   virtual void do_symbolic_factorization(CSCSymmetricMatrix& matrix) = 0;
   virtual void do_numerical_factorization(CSCSymmetricMatrix& matrix) = 0;
   virtual std::vector<double> solve(CSCSymmetricMatrix& matrix, const std::vector<double>& rhs) = 0;

   [[nodiscard]] virtual std::tuple<int, int, int> get_inertia() const = 0;
   [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
   [[nodiscard]] virtual bool matrix_is_singular() const = 0;
   [[nodiscard]] virtual int rank() const = 0;
};

template <class MatrixType>
class LinearSolverTest {

};

#endif // LINEARSOLVER_H
