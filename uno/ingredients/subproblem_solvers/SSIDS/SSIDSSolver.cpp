// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "SSIDSSolver.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   SSIDSSolver::SSIDSSolver(): DirectSymmetricIndefiniteLinearSolver<double>() {
      spral_ssids_default_options(&this->workspace.options);
      this->workspace.options.array_base = 1; // Fortran indexing
      this->workspace.options.print_level = -1; // no printing
   }

   void SSIDSSolver::initialize_memory() {
      this->workspace.n = static_cast<int>(this->linear_system.dimension);
      this->workspace.nnz = static_cast<int>(this->linear_system.number_nonzeros);
   }

   void SSIDSSolver::do_symbolic_analysis() {
      assert(!this->analysis_performed);

      spral_ssids_analyse_coord(this->workspace.n, nullptr, this->workspace.nnz, this->linear_system.matrix_row_indices.data(),
         this->linear_system.matrix_column_indices.data(), nullptr, &this->workspace.akeep, &this->workspace.options,
         &this->workspace.inform);
      if (this->workspace.inform.flag < 0) {
         spral_ssids_free(&this->workspace.akeep, &this->workspace.fkeep);
         throw std::runtime_error("SSIDS could not compute the symbolic analysis");
      }
      this->analysis_performed = true;
   }

   void SSIDSSolver::do_numerical_factorization(bool is_matrix_positive_definite) {
      assert(this->analysis_performed);

      spral_ssids_factor(is_matrix_positive_definite, nullptr, nullptr, this->linear_system.matrix_values.data(), nullptr,
         this->workspace.akeep, &this->workspace.fkeep, &this->workspace.options, &this->workspace.inform);
      if(this->workspace.inform.flag < 0) {
         spral_ssids_free(&this->workspace.akeep, &this->workspace.fkeep);
         throw std::runtime_error("SSIDS could not compute the factorization");
      }
      this->factorization_performed = true;
   }

   void SSIDSSolver::solve_indefinite_system(double* solution) {
      assert(this->factorization_performed);

      // copy rhs into solution (overwritten by SSIDS)
      const size_t dimension = static_cast<size_t>(this->workspace.n);
      view(solution, dimension) = this->linear_system.rhs.view();

      spral_ssids_solve1(0, solution, this->workspace.akeep, this->workspace.fkeep, &this->workspace.options,
         &this->workspace.inform);
      if (this->workspace.inform.flag < 0) {
         spral_ssids_free(&this->workspace.akeep, &this->workspace.fkeep);
         throw std::runtime_error("SSIDS could not solve the linear system");
      }
   }

   void SSIDSSolver::solve_indefinite_system(const double* rhs, double* solution, size_t number_of_rhs) {
      assert(this->factorization_performed);

      // copy the rhs block into the solution block (overwritten by SSIDS); both are column-major with leading dimension n
      const size_t dimension = static_cast<size_t>(this->workspace.n);
      view(solution, number_of_rhs * dimension) = view(rhs, number_of_rhs * dimension);

      spral_ssids_solve(0, static_cast<int>(number_of_rhs), solution, this->workspace.n, this->workspace.akeep,
         this->workspace.fkeep, &this->workspace.options, &this->workspace.inform);
      if (this->workspace.inform.flag < 0) {
         spral_ssids_free(&this->workspace.akeep, &this->workspace.fkeep);
         throw std::runtime_error("SSIDS could not solve the linear system");
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

   LinearSystem& SSIDSSolver::get_linear_system() {
      return this->linear_system;
   }

   COOLinearSystem& SSIDSSolver::get_coo_linear_system() {
      return this->linear_system;
   }
} // namespace