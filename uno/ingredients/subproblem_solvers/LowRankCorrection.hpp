// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LOWRANKCORRECTION_H
#define UNO_LOWRANKCORRECTION_H

#include <optional>
#include <vector>
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
   // - solve(r, x): solves the system with the full matrix, given the factorization of A

   // the augmented matrix contains the full Hessian: plain solve
   struct NoLowRankCorrection {
      static constexpr bool uses_explicit_hessian_matrix = true;
      void update(const Subproblem& /*subproblem*/, DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/) { }
      void solve(DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, const Vector<double>& rhs,
         Vector<double>& solution);
   };

   // The Hessian approximation is a low-rank correction to a diagonal part: H = δ I + E C Eᵀ, C = diag(c_i).
   // The augmented matrix
   // (δ I + E C Eᵀ     Jᵀ) = (δ I   Jᵀ) + (E) C (Eᵀ  0)    (*)
   // (J                0 )   (J     0 )   (0)
   // is inverted with the Woodbury formula: A contains the diagonal part only, and
   // (*)⁻¹ r = A⁻¹ r - A⁻¹ E T⁻¹ Eᵀ A⁻¹ r,  with T = C⁻¹ + Eᵀ A⁻¹ E.
   // With the split A⁻¹ = Rᵀ M R exposed by the linear solver (R = (P L Pᵀ)⁻¹, M = (P D Pᵀ)⁻¹), the symmetric form
   // W = R E, Z = M W, T = C⁻¹ + Wᵀ Z only requires forward and diagonal solves per correction column (no backward
   // solve), and (*)⁻¹ r = Rᵀ M (u - W T⁻¹ Zᵀ u) with u = R r costs a single full solve per right-hand side.
   // Without partial solves (R = I, M = A⁻¹), this is the standard Woodbury formula.
   // W, Z and the Bunch-Kaufman factors of T are computed once per factorization of A and reused for every
   // right-hand side (direction, SOC)
   class WoodburyCorrection {
   public:
      static constexpr bool uses_explicit_hessian_matrix = false;

      explicit WoodburyCorrection(const DirectQuasiNewtonHessian& hessian_model);
      void update(const Subproblem& subproblem, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver);
      void solve(DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, const Vector<double>& rhs,
         Vector<double>& solution);

   protected:
      const DirectQuasiNewtonHessian& hessian_model;
      size_t correction_rank{0};
      // buffers, reallocated only when their dimensions change
      std::optional<DenseMatrix<double>> W{}; // R E
      std::optional<DenseMatrix<double>> Z{}; // M R E
      std::optional<DenseMatrix<double>> T{}; // Bunch-Kaufman factors of C⁻¹ + Wᵀ Z
      std::vector<int> ipiv{};
      Vector<double> c{};
      Vector<double> d{};
   };
} // namespace

#endif // UNO_LOWRANKCORRECTION_H