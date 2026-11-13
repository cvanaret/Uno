// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
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
   void NoLowRankCorrection::solve(DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, const Vector<double>& rhs,
         Vector<double>& solution) {
      linear_solver.solve_indefinite_system(rhs.data(), solution.data(), 1);
   }

   namespace {
      // (re)allocate a matrix only when its dimensions change
      void resize(std::optional<DenseMatrix<double>>& matrix, size_t number_rows, size_t number_columns) {
         if (!matrix.has_value() || matrix->number_rows != number_rows || matrix->number_columns != number_columns) {
            matrix.emplace(number_rows, number_columns);
         }
      }
   } // namespace

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
      const size_t dimension = linear_solver.get_linear_system().dimension;
      resize(this->W, dimension, this->correction_rank);
      resize(this->Z, dimension, this->correction_rank);
      resize(this->T, this->correction_rank, this->correction_rank);
      this->c.resize(this->correction_rank);
      this->d.resize(this->correction_rank);

      // W = E: the correction columns live in the primal rows, the other rows are 0
      const size_t number_variables = subproblem.problem.model.number_variables;
      for (size_t column_index: Range(this->correction_rank)) {
         auto W_column = this->W->column(column_index);
         view(W_column, 0, number_variables) = this->hessian_model.get_correction_column(column_index);
         view(W_column, number_variables, dimension).fill(0.);
      }
      // W = R E, then Z = M W (all correction columns as right-hand sides at once)
      linear_solver.solve_lower_factor(this->W->data(), this->correction_rank);
      std::copy(this->W->data(), this->W->data() + dimension * this->correction_rank, this->Z->data());
      linear_solver.solve_middle_factor(this->Z->data(), this->correction_rank);

      // T = C⁻¹ + Wᵀ Z, then factorize it in place
      this->T->fill(0.);
      for (size_t column_index: Range(this->correction_rank)) {
         this->T->entry(column_index, column_index) = 1./this->hessian_model.get_correction_column_scaling(column_index);
      }
      *this->T += transpose(*this->W) * *this->Z;
      DEBUG2 << "T = " << *this->T;
      auto [success, ipiv] = this->T->compute_bunch_kaufman_factorization();
      if (!success) {
         throw std::runtime_error("WoodburyCorrection: the Bunch-Kaufman factorization failed");
      }
      this->ipiv = std::move(ipiv);
   }

   // solution = Rᵀ M (u - W T⁻¹ Zᵀ u) with u = R rhs
   void WoodburyCorrection::solve(DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, const Vector<double>& rhs,
         Vector<double>& solution) {
      if (this->correction_rank == 0) {
         linear_solver.solve_indefinite_system(rhs.data(), solution.data(), 1);
         return;
      }
      solution = rhs;
      linear_solver.solve_lower_factor(solution.data(), 1);
      this->c = transpose(*this->Z) * solution;
      if (!solve_bunch_kaufman(*this->T, this->c, this->d, this->ipiv)) {
         throw std::runtime_error("WoodburyCorrection: the Bunch-Kaufman solve failed");
      }
      solution -= *this->W * this->d;
      linear_solver.solve_middle_factor(solution.data(), 1);
      linear_solver.solve_upper_factor(solution.data(), 1);
   }
} // namespace