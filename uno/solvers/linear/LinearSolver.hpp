#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <vector>
#include <linear_algebra/SymmetricMatrix.hpp>

class LinearSolver {
public:
   explicit LinearSolver(size_t max_dimension): max_dimension(max_dimension) {};
   virtual ~LinearSolver() = default;

   // matrix is not declared const, since Fortran-based solvers may need to temporarily reindex the coordinates
   virtual void factorize(size_t dimension, const SymmetricMatrix& matrix) = 0;
   virtual void do_symbolic_factorization(size_t dimension, const SymmetricMatrix& matrix) = 0;
   virtual void do_numerical_factorization(size_t dimension, const SymmetricMatrix& matrix) = 0;
   virtual void solve(size_t dimension, const SymmetricMatrix& matrix, const std::vector<double>& rhs, std::vector<double>& result) = 0;

   [[nodiscard]] virtual std::tuple<size_t, size_t, size_t> get_inertia() const = 0;
   [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
   [[nodiscard]] virtual bool matrix_is_singular() const = 0;
   [[nodiscard]] virtual size_t rank() const = 0;

protected:
   const size_t max_dimension;
};

#endif // LINEARSOLVER_H
