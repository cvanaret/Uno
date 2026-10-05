// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <stdexcept>
#include "LowRankCorrection.hpp"
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/hessian_models/quasi_newton/direct/DirectQuasiNewtonHessian.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "model/Model.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Transpose.hpp"
#include "tools/Logger.hpp"

namespace uno {
   WoodburyCorrection::WoodburyCorrection(const DirectQuasiNewtonHessian& hessian_model): hessian_model(hessian_model) {
      assert(!this->hessian_model.has_hessian_matrix());
   }

   // precondition: the augmented matrix A has been factorized and is nonsingular
   void WoodburyCorrection::update(const Subproblem& subproblem, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver) {
      this->correction_rank = this->hessian_model.get_correction_rank();
      DEBUG << "Correction rank: " << this->correction_rank << '\n';
      if (this->correction_rank == 0) {
         return;
      }
      const size_t dimension = subproblem.number_variables + subproblem.number_constraints;
      // assemble the correction columns into E (E is taller than the correction matrix; the extra rows stay 0)
      this->E.emplace(dimension, this->correction_rank);
      this->H.emplace(dimension, this->correction_rank);
      for (size_t column_index: Range(this->correction_rank)) {
         const auto correction_column = this->hessian_model.get_correction_column(column_index);
         for (size_t row_index: Range(subproblem.problem.model.number_variables)) {
            this->E->entry(row_index, column_index) = correction_column[row_index];
         }
      }
      // solve A H = E with all correction columns as right-hand sides at once (column-major blocks)
      linear_solver.solve_indefinite_system(this->E->data(), this->H->data(), this->correction_rank);
      // T = P⁻¹ + Eᵀ H, then factorize it in place
      this->T.emplace(this->correction_rank, this->correction_rank);
      for (size_t column_index: Range(this->correction_rank)) {
         this->T->entry(column_index, column_index) = 1./this->hessian_model.get_correction_column_scaling(column_index);
      }
      *this->T += transpose(*this->E) * *this->H;
      DEBUG2 << "T = " << *this->T;
      auto [success, ipiv] = this->T->compute_bunch_kaufman_factorization();
      if (!success) {
         throw std::runtime_error("WoodburyCorrection: the Bunch-Kaufman factorization failed");
      }
      this->ipiv = std::move(ipiv);
   }

   // b := b - H T⁻¹ Eᵀ b
   void WoodburyCorrection::apply(Vector<double>& b) const {
      if (this->correction_rank == 0) {
         return;
      }
      Vector<double> c(this->correction_rank);
      Vector<double> d(this->correction_rank);
      c = transpose(*this->E) * b;
      if (!solve_bunch_kaufman(*this->T, c, d, this->ipiv)) {
         throw std::runtime_error("WoodburyCorrection: the Bunch-Kaufman solve failed");
      }
      b -= *this->H * d;
   }
} // namespace