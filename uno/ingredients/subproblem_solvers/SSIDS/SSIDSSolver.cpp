// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "SSIDSSolver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"

namespace uno {
   SSIDSSolver::SSIDSSolver(): DirectSymmetricIndefiniteLinearSolver<double>() {
      spral_ssids_default_options(&this->workspace.options);
      this->workspace.options.array_base = 1; // Fortran indexing
      this->workspace.options.print_level = -1; // no printing
   }

   void SSIDSSolver::initialize_hessian(const Subproblem& subproblem) {
      this->coo_workspace.initialize_hessian(subproblem);
      this->workspace.n = static_cast<int>(subproblem.number_variables);
      this->workspace.nnz = static_cast<int>(this->coo_workspace.number_matrix_nonzeros);
   }

   void SSIDSSolver::initialize_augmented_system(const Subproblem& subproblem) {
      this->coo_workspace.initialize_augmented_system(subproblem);
      this->workspace.n = static_cast<int>(subproblem.number_variables + subproblem.number_constraints);
      this->workspace.nnz = static_cast<int>(this->coo_workspace.number_matrix_nonzeros);
   }

   void SSIDSSolver::do_symbolic_analysis() {
      spral_ssids_analyse_coord(this->workspace.n, nullptr, this->workspace.nnz, this->coo_workspace.matrix_row_indices.data(),
         this->coo_workspace.matrix_column_indices.data(), nullptr, &this->workspace.akeep, &this->workspace.options,
         &this->workspace.inform);
      if (this->workspace.inform.flag < 0) {
         spral_ssids_free(&this->workspace.akeep, &this->workspace.fkeep);
         throw std::runtime_error("SSIDS could not compute the symbolic analysis");
      }
   }

   void SSIDSSolver::do_numerical_factorization(const double* matrix_values, bool is_matrix_positive_definite) {
      spral_ssids_factor(is_matrix_positive_definite, nullptr, nullptr, matrix_values, nullptr, this->workspace.akeep,
         &this->workspace.fkeep, &this->workspace.options, &this->workspace.inform);
      if(this->workspace.inform.flag < 0) {
         spral_ssids_free(&this->workspace.akeep, &this->workspace.fkeep);
         throw std::runtime_error("SSIDS could not compute the factorization");
      }
   }

   void SSIDSSolver::solve_indefinite_system(const Vector<double>& /*matrix_values*/, const Vector<double>& rhs, Vector<double>& result) {
      result = rhs;
      spral_ssids_solve1(0, result.data(), this->workspace.akeep, this->workspace.fkeep, &this->workspace.options,
         &this->workspace.inform);
      if (this->workspace.inform.flag < 0) {
         spral_ssids_free(&this->workspace.akeep, &this->workspace.fkeep);
         throw std::runtime_error("SSIDS could not solve the linear system");
      }
   }

   void SSIDSSolver::solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) {
      // set up the linear system by evaluating the functions at the current iterate
      this->coo_workspace.set_up_linear_system(statistics, subproblem, *this, current_evaluations, warmstart_information);
      // solve the linear system
      this->solve_indefinite_system(this->coo_workspace.matrix_values, this->coo_workspace.rhs, this->coo_workspace.solution);
      // assemble the full primal-dual direction
      subproblem.assemble_primal_dual_direction(this->coo_workspace.solution, direction);
      if (this->matrix_is_singular()) {
         direction.status = SubproblemStatus::INFEASIBLE;
      }
   }

   Inertia SSIDSSolver::get_inertia() const {
      // rank = number_positive_eigenvalues + number_negative_eigenvalues
      // n = rank + number_zero_eigenvalues
      const size_t rank = this->rank();
      const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
      const size_t number_positive_eigenvalues = rank - number_negative_eigenvalues;
      const size_t number_zero_eigenvalues = static_cast<size_t>(this->workspace.n) - rank;
      return {number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues};
   }

   size_t SSIDSSolver::number_negative_eigenvalues() const {
      return static_cast<size_t>(this->workspace.inform.num_neg);
   }

   bool SSIDSSolver::matrix_is_singular() const {
      return (this->workspace.inform.matrix_rank < this->workspace.n);
   }

   size_t SSIDSSolver::rank() const {
      return static_cast<size_t>(this->workspace.inform.matrix_rank);
   }

   SolverWorkspace& SSIDSSolver::get_workspace() {
      return this->coo_workspace;
   }
} // namespace