// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <cmath>
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
   namespace {
      Inertia bunch_kaufman_inertia(const DenseMatrix<double>& factors, const std::vector<int>& ipiv) {
         size_t positive = 0, negative = 0, zero = 0;
         const auto count = [&](double value) {
            if (!std::isfinite(value)) throw std::runtime_error("WoodburyCorrection: non-finite pivot");
            if (0. < value) ++positive;
            else if (value < 0.) ++negative;
            else ++zero;
         };
         for (size_t k = 0; k < ipiv.size();) {
            if (0 < ipiv[k]) {
               count(factors.entry(k, k));
               ++k;
            }
            else {
               if (ipiv[k] == 0 || k + 1 == ipiv.size() || ipiv[k + 1] != ipiv[k]) {
                  throw std::runtime_error("WoodburyCorrection: invalid Bunch-Kaufman pivots");
               }
               const double a = factors.entry(k, k), b = factors.entry(k + 1, k), c = factors.entry(k + 1, k + 1);
               const double scale = std::max({std::abs(a), std::abs(b), std::abs(c)});
               if (scale == 0.) {
                  zero += 2;
               }
               else {
                  const double x = a / scale, y = b / scale, z = c / scale;
                  const double center = (x + z) / 2.;
                  const double large = center + std::copysign(std::hypot((x - z) / 2., y), center);
                  count(large);
                  count(large == 0. ? 0. : std::fma(x, z, -y * y) / large);
               }
               k += 2;
            }
         }
         return {positive, negative, zero};
      }
   } // namespace
   WoodburyCorrection::WoodburyCorrection(const DirectQuasiNewtonHessian& hessian_model): hessian_model(hessian_model) {
      assert(!this->hessian_model.has_hessian_matrix());
   }

   // precondition: the augmented matrix A has been factorized and is nonsingular
   void WoodburyCorrection::update(const Subproblem& subproblem, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver) {
      this->correction_rank = this->hessian_model.get_correction_rank();
      this->correction_inertia.reset();
      this->scaling_positive = this->scaling_negative = 0;
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
      // with the correction E C Eᵀ (C = diag(c_i), c_i the coefficients of the rank-one terms): T = C⁻¹ + Eᵀ H, then
      // factorize it in place
      this->T.emplace(this->correction_rank, this->correction_rank);
      for (size_t column_index: Range(this->correction_rank)) {
         const double coefficient = this->hessian_model.get_correction_column_scaling(column_index);
         const double inverse = 1. / coefficient;
         if (!std::isfinite(coefficient) || !std::isfinite(inverse)) {
            throw std::runtime_error("WoodburyCorrection: invalid rank-one coefficient");
         }
         if (0. < coefficient) ++this->scaling_positive;
         else ++this->scaling_negative;
         this->T->entry(column_index, column_index) = inverse;
      }
      *this->T += transpose(*this->E) * *this->H;
      for (size_t column: Range(this->correction_rank)) {
         for (size_t row: Range(column, this->correction_rank)) {
            if (!std::isfinite(this->T->entry(row, column))) {
               throw std::runtime_error("WoodburyCorrection: non-finite correction matrix");
            }
         }
      }
      DEBUG2 << "T = " << *this->T;
      auto [success, ipiv] = this->T->compute_bunch_kaufman_factorization();
      this->correction_inertia.emplace(bunch_kaufman_inertia(*this->T, ipiv));
      if (!success && this->correction_inertia->zero == 0) {
         throw std::runtime_error("WoodburyCorrection: the Bunch-Kaufman factorization failed");
      }
      this->ipiv = std::move(ipiv);
   }

   Inertia WoodburyCorrection::get_inertia(const Inertia& base_inertia) const {
      if (this->correction_rank == 0) return base_inertia;
      assert(base_inertia.zero == 0 && this->correction_inertia);
      // Haynsworth: inertia(A + E C Eᵀ) = inertia(A) + inertia(-T) - inertia(-C⁻¹).
      const size_t positive = base_inertia.positive + this->correction_inertia->negative;
      const size_t negative = base_inertia.negative + this->correction_inertia->positive;
      if (positive < this->scaling_negative || negative < this->scaling_positive) {
         throw std::runtime_error("WoodburyCorrection: inconsistent inertia");
      }
      return {positive - this->scaling_negative, negative - this->scaling_positive, this->correction_inertia->zero};
   }

   void WoodburyCorrection::add_matrix_vector_product(View<const double> vector, View<double> result) const {
      for (size_t column_index: Range(this->correction_rank)) {
         const auto column = this->hessian_model.get_correction_column(column_index);
         const double coefficient = this->hessian_model.get_correction_column_scaling(column_index) *
            dot(column, view(vector.data(), column.size()));
         blas1::add(column.size(), coefficient, column.data(), result.data());
      }
   }

   // b := b - H T⁻¹ Eᵀ b
   void WoodburyCorrection::apply(Vector<double>& b) const {
      if (this->correction_rank == 0) {
         return;
      }
      if (!this->correction_inertia || 0 < this->correction_inertia->zero) {
         throw std::runtime_error("WoodburyCorrection: the full matrix is singular");
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