// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LOWRANKCORRECTION_H
#define UNO_LOWRANKCORRECTION_H

#include <optional>
#include <vector>
#include "ingredients/inertia_correction_strategies/Inertia.hpp"
#include "linear_algebra/DenseMatrix.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class DirectQuasiNewtonHessian;
   class Subproblem;

   // Policies for EQPSolver. Interface:
   // - uses_explicit_hessian_matrix: whether the augmented matrix contains the full Hessian
   // - update(): called once per numerical factorization of the augmented matrix A
   // - get_inertia(): returns the inertia including the low-rank terms
   // - apply(b): turns b = A⁻¹ r into the solution of the system with the full matrix

   // the augmented matrix contains the full Hessian: nothing to do
   struct NoLowRankCorrection {
      static constexpr bool uses_explicit_hessian_matrix = true;
      void update(const Subproblem& /*subproblem*/, DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/) { }
      void apply(Vector<double>& /*b*/) const { }
      [[nodiscard]] Inertia get_inertia(const Inertia& base_inertia) const { return base_inertia; }
   };

   // The Hessian approximation is a low-rank correction to a diagonal part: H = δ I + E P Eᵀ.
   // The augmented matrix
   // (δ I + E P Eᵀ     Jᵀ) = (δ I   Jᵀ) + (E) P (Eᵀ  0)    (*)
   // (J                0 )   (J     0 )   (0)
   // is inverted with the Woodbury formula: A contains the diagonal part only, and
   // (*)⁻¹ r = b - H (P⁻¹ + Eᵀ H)⁻¹ Eᵀ b,  with b = A⁻¹ r and H = A⁻¹ E
   // H and the Bunch-Kaufman factors of T = P⁻¹ + Eᵀ H are computed once per factorization of A and reused
   // for every right-hand side (direction, SOC)
   class WoodburyCorrection {
   public:
      static constexpr bool uses_explicit_hessian_matrix = false;

      explicit WoodburyCorrection(const DirectQuasiNewtonHessian& hessian_model);
      void update(const Subproblem& subproblem, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver);
      void apply(Vector<double>& b) const;
      [[nodiscard]] Inertia get_inertia(const Inertia& base_inertia) const;

   protected:
      const DirectQuasiNewtonHessian& hessian_model;
      size_t correction_rank{0};
      std::optional<DenseMatrix<double>> E{};
      std::optional<DenseMatrix<double>> H{};
      std::optional<DenseMatrix<double>> T{}; // Bunch-Kaufman factors of P⁻¹ + Eᵀ H
      std::vector<int> ipiv{};
      std::optional<Inertia> correction_inertia{};
      size_t scaling_positive{0}, scaling_negative{0};
   };
} // namespace

#endif // UNO_LOWRANKCORRECTION_H