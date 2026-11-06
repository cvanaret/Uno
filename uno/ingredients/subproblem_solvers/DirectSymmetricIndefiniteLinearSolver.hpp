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

      [[nodiscard]] virtual Inertia get_inertia() const = 0;
      [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
      // [[nodiscard]] virtual bool matrix_is_positive_definite() const = 0;
      [[nodiscard]] virtual bool matrix_is_singular() const = 0;
      [[nodiscard]] virtual size_t rank() const = 0;

      // increase the pivoting tolerance (to improve the accuracy of the next factorizations) if the solver supports it.
      // Returns false if the tolerance cannot be increased (not supported or maximum value reached).
      // The new tolerance is used from the next numerical factorization on, and is never decreased
      [[nodiscard]] virtual bool increase_pivoting_tolerance() { return false; }

      // split solves with the factors of A = P L D Lᵀ Pᵀ: A⁻¹ = (P Lᵀ Pᵀ)⁻¹ (P D Pᵀ)⁻¹ (P L Pᵀ)⁻¹. Each function overwrites
      // `number_of_rhs` column-major right-hand sides. Solvers without partial solves treat A as the middle factor
      // (L = I), so that the composition is always the full solve
      virtual void solve_lower_factor(ElementType* /*rhs*/, size_t /*number_of_rhs*/) { } // (P L Pᵀ)⁻¹
      virtual void solve_middle_factor(ElementType* rhs, size_t number_of_rhs) { // (P D Pᵀ)⁻¹
         this->solve_indefinite_system(rhs, rhs, number_of_rhs);
      }
      virtual void solve_upper_factor(ElementType* /*rhs*/, size_t /*number_of_rhs*/) { } // (P Lᵀ Pᵀ)⁻¹
   };
} // namespace

#endif // UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H