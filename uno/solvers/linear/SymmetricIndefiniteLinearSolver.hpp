// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_SYMMETRICINDEFINITELINEARSOLVER_H

#include <vector>
#include <linear_algebra/SymmetricMatrix.hpp>

template <typename T>
class SymmetricIndefiniteLinearSolver {
public:
   explicit SymmetricIndefiniteLinearSolver(size_t max_dimension): max_dimension(max_dimension) {};
   virtual ~SymmetricIndefiniteLinearSolver() = default;

   virtual void factorize(const SymmetricMatrix<T>& matrix) = 0;
   virtual void do_symbolic_factorization(const SymmetricMatrix<T>& matrix) = 0;
   virtual void do_numerical_factorization(const SymmetricMatrix<T>& matrix) = 0;
   virtual void solve_indefinite_system(const SymmetricMatrix<T>& matrix, const std::vector<T>& rhs, std::vector<T>& result) = 0;

   [[nodiscard]] virtual std::tuple<size_t, size_t, size_t> get_inertia() const = 0;
   [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
   // [[nodiscard]] virtual bool matrix_is_positive_definite() const = 0;
   [[nodiscard]] virtual bool matrix_is_singular() const = 0;
   [[nodiscard]] virtual size_t rank() const = 0;

protected:
   const size_t max_dimension;
};

#endif // UNO_SYMMETRICINDEFINITELINEARSOLVER_H
