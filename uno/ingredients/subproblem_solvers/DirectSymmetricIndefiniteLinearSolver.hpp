// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H

#include <cstddef>
#include "SymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/inertia_correction_strategies/Inertia.hpp"

namespace uno {
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver: public SymmetricIndefiniteLinearSolver<ElementType> {
   public:
      DirectSymmetricIndefiniteLinearSolver() = default;
      ~DirectSymmetricIndefiniteLinearSolver() override = default;

      virtual void do_symbolic_analysis() = 0;
      virtual void do_numerical_factorization(bool is_matrix_positive_definite) = 0;

      // hint, set before the symbolic analysis, for solvers that need the expected sign of each pivot
      // (e.g. HiPO): the leading expected_inertia.positive pivots are positive, the next
      // expected_inertia.negative ones are negative. No-op for solvers that do not need it.
      virtual void set_expected_inertia(const Inertia& /*expected_inertia*/) { }

      [[nodiscard]] virtual Inertia get_inertia() const = 0;
      [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
      // [[nodiscard]] virtual bool matrix_is_positive_definite() const = 0;
      [[nodiscard]] virtual bool matrix_is_singular() const = 0;
      [[nodiscard]] virtual size_t rank() const = 0;
   };
} // namespace

#endif // UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H