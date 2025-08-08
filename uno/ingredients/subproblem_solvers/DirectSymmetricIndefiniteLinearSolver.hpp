// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H

#include "SymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/regularization_strategies/Inertia.hpp"

namespace uno {
   template <typename IndexType, typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver: public SymmetricIndefiniteLinearSolver<IndexType, ElementType> {
   public:
      DirectSymmetricIndefiniteLinearSolver() = default;
      ~DirectSymmetricIndefiniteLinearSolver() override = default;

      virtual void do_symbolic_analysis() = 0;
      virtual void do_numerical_factorization(const double* matrix_values) = 0;

      [[nodiscard]] virtual Inertia get_inertia() const = 0;
      [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
      // [[nodiscard]] virtual bool matrix_is_positive_definite() const = 0;
      [[nodiscard]] virtual bool matrix_is_singular() const = 0;
      [[nodiscard]] virtual size_t rank() const = 0;
   };
} // namespace

#endif // UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H