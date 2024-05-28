// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DIRECTINDEFINITELINEARSOLVER_H
#define UNO_DIRECTINDEFINITELINEARSOLVER_H

#include "solvers/linear/SymmetricIndefiniteLinearSolver.hpp"

template <typename IndexType, typename ElementType>
class DirectIndefiniteLinearSolver: public SymmetricIndefiniteLinearSolver<IndexType, ElementType> {
public:
   explicit DirectIndefiniteLinearSolver(size_t max_dimension): SymmetricIndefiniteLinearSolver<IndexType, ElementType>(max_dimension) {
   }
   virtual ~DirectIndefiniteLinearSolver() = default;

   virtual void do_symbolic_factorization(const SymmetricMatrix<IndexType, ElementType>& matrix) = 0;
   virtual void do_numerical_factorization(const SymmetricMatrix<IndexType, ElementType>& matrix) = 0;
   
   [[nodiscard]] virtual std::tuple<size_t, size_t, size_t> get_inertia() const = 0;
   [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
   // [[nodiscard]] virtual bool matrix_is_positive_definite() const = 0;
   [[nodiscard]] virtual bool matrix_is_singular() const = 0;
   [[nodiscard]] virtual size_t rank() const = 0;
};

#endif // UNO_DIRECTINDEFINITELINEARSOLVER_H