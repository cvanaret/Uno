// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "WoodburyEQPSolver.hpp"
#include "LinearSystem.hpp"
#include "SymmetricIndefiniteLinearSolverFactory.hpp"
#include "ingredients/hessian_models/quasi_newton/direct/DirectQuasiNewtonHessian.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Transpose.hpp"

namespace uno {
   WoodburyEQPSolver::WoodburyEQPSolver(const DirectQuasiNewtonHessian& hessian_model, const Options& options):
         SubproblemSolver(),
         hessian_model(hessian_model),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"))) {
      assert(!this->hessian_model.has_hessian_matrix());
   }

   void WoodburyEQPSolver::initialize_memory(const Subproblem& subproblem) {
      // access the linear system of the linear solver
      auto& linear_system = this->linear_solver->get_linear_system();
      linear_system.initialize_augmented_system(subproblem);
      this->linear_solver->initialize_memory();
   }

   void WoodburyEQPSolver::solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& /*initial_point*/, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      if (is_finite(trust_region_radius)) {
         throw std::runtime_error("The direct linear solver does not support a trust region");
      }

      // access the linear system
      auto& linear_system = this->linear_solver->get_linear_system();

      // set up the linear system by evaluating the functions at the current iterate
      if (warmstart_information.new_iterate) {
         // assemble the augmented matrix with only the diagonal part of the quasi-Newton Hessian
         subproblem.evaluate_lagrangian_hessian(statistics, linear_system.matrix_values.data());
         const size_t number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
         subproblem.evaluate_jacobian(linear_system.matrix_values.data() + number_hessian_nonzeros, current_evaluations);

         // perform the symbolic analysis once and for all
         if (!this->analysis_performed) {
            DEBUG << "Performing symbolic analysis of the indefinite system\n";
            this->linear_solver->do_symbolic_analysis();
            this->analysis_performed = true;
         }

         // regularize the augmented matrix (this calls the analysis and the factorization)
         subproblem.regularize_augmented_matrix(statistics, linear_system.matrix_values.data(),
            subproblem.dual_regularization_factor(), *this->linear_solver);

         // assemble the RHS
         subproblem.assemble_augmented_rhs(current_evaluations, linear_system.rhs);
      }

      // solve the linear system with only the diagonal part and store the result in solution_diagonal_part
      Vector<double> solution_diagonal_part(subproblem.number_variables + subproblem.number_constraints);
      this->linear_solver->solve_indefinite_system(solution_diagonal_part.data());
      if (this->linear_solver->matrix_is_singular()) {
         direction.status = SubproblemStatus::INFEASIBLE;
         return;
      }

      // compute the low-rank correction
      this->compute_low_rank_correction(subproblem, linear_system, solution_diagonal_part);

      // assemble the full primal-dual direction
      subproblem.assemble_primal_dual_direction(solution_diagonal_part, direction);
   }

   SolverWorkspace& WoodburyEQPSolver::get_workspace() {
      return this->linear_solver->get_linear_system();
   }

   // protected members

   void WoodburyEQPSolver::compute_low_rank_correction(const Subproblem& subproblem, LinearSystem& linear_system, Vector<double>& b) const {
      DEBUG << "b = " << b << '\n';
      const size_t correction_rank = this->hessian_model.get_correction_rank();
      DEBUG << "Correction rank: " << correction_rank << '\n';
      if (0 < correction_rank) {
         // compute correction_rank backsolves with the correction columns as RHS
         DenseMatrix<double> E(subproblem.number_variables + subproblem.number_constraints, correction_rank);
         DenseMatrix<double> H(subproblem.number_variables + subproblem.number_constraints, correction_rank);
         for (size_t column_index: Range(correction_rank)) {
            const auto correction_column = this->hessian_model.get_correction_column(column_index);
            // copy each correction column into E (note: E is higher than the correction matrix)
            for (size_t row_index: Range(subproblem.problem.model.number_variables)) {
               E.entry(row_index, column_index) = correction_column[row_index];
            }
            // solve a linear system A H_j = E_j
            linear_system.rhs = E.column(column_index);
            this->linear_solver->solve_indefinite_system(H.column(column_index).data());
         }
         DEBUG2 << "E = " << E;
         DEBUG2 << "H = " << H;
         // compute c = Eᵀ b
         Vector<double> c(correction_rank);
         c = transpose(E)*b; // TODO move to constructor
         DEBUG2 << "c = " << c << '\n';
         // construct T = P⁻¹ + Eᵀ H: first set P⁻¹ into T, then add Eᵀ H
         DenseMatrix<double> T(correction_rank, correction_rank);
         for (size_t column_index: Range(correction_rank)) {
            T.entry(column_index, column_index) = 1./this->hessian_model.get_correction_column_scaling(column_index);
         }
         T += transpose(E)*H;
         DEBUG2 << "T = " << T;
         // solve T d = c by computing a Bunch-Kaufman factorization of T
         Vector<double> d(correction_rank);
         const bool success = WoodburyEQPSolver::solve_dense_indefinite_system(T, c, d);
         DEBUG2 << "Bunch-Kaufman success: " << success << '\n';
         DEBUG2 << "d = " << d << '\n';
         // add the correction to b: b := b - H d
         b -= H * d;
         DEBUG2 << "b = " << b << '\n';
      }
   }

   bool WoodburyEQPSolver::solve_dense_indefinite_system(DenseMatrix<double>& T, const Vector<double>& c, Vector<double>& d) {
      auto [success, ipiv] = T.compute_bunch_kaufman_factorization();
      if (success) {
         success = solve_bunch_kaufman(T, c, d, ipiv);
      }
      return success;
   }
} // namespace