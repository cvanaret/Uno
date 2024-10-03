// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H

#include "solvers/SymmetricIndefiniteLinearSolver.hpp"

namespace uno {
   template <typename IndexType, typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver: public SymmetricIndefiniteLinearSolver<IndexType, ElementType> {
   public:
      explicit DirectSymmetricIndefiniteLinearSolver(size_t dimension) : SymmetricIndefiniteLinearSolver<IndexType, ElementType>(dimension) { };
      virtual ~DirectSymmetricIndefiniteLinearSolver() = default;

      virtual void factorize(const SymmetricMatrix<IndexType, ElementType>& matrix) = 0;
      virtual void do_symbolic_factorization(const SymmetricMatrix<IndexType, ElementType>& matrix) = 0;
      virtual void do_numerical_factorization(const SymmetricMatrix<IndexType, ElementType>& matrix) = 0;

      [[nodiscard]] virtual std::tuple<size_t, size_t, size_t> get_inertia() const = 0;
      [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
      // [[nodiscard]] virtual bool matrix_is_positive_definite() const = 0;
      [[nodiscard]] virtual bool matrix_is_singular() const = 0;
      [[nodiscard]] virtual size_t rank() const = 0;
   };
} // namespace

#endif // UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H