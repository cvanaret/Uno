// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "LBFGSHessian.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Range.hpp"

#ifdef HAS_LAPACK
#include "fortran_interface.h"
#define LAPACK_cholesky_factorization FC_GLOBAL_(dpotrf, DPOTRF)
extern "C" {
   void LAPACK_cholesky_factorization(char* uplo, int* n, double* a, int* lda, int* info);
}
#endif

namespace uno {
   LBFGSHessian::LBFGSHessian(size_t dimension, size_t memory_size):
         HessianModel(),
         dimension(dimension),
         memory_size(memory_size),
         S_matrix(dimension, memory_size),
         Y_matrix(dimension, memory_size),
         L_matrix(memory_size, memory_size),
         D_matrix(memory_size),
         M_matrix(memory_size, memory_size) {
   }

   void LBFGSHessian::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) const {
      // do nothing
   }

   void LBFGSHessian::notify_accepted_iterate(const Iterate& current_iterate, const Iterate& trial_iterate) {
      this->update_memory(current_iterate, trial_iterate);
      this->hessian_recomputation_required = true;
   }

   void LBFGSHessian::evaluate_hessian(Statistics& /*statistics*/, const Model& /*model*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, SymmetricMatrix<size_t, double>& /*hessian*/) {
      if (this->hessian_recomputation_required) {
         this->compute_hessian_representation();
         this->hessian_recomputation_required = false;
      }

      // TODO
      throw std::runtime_error("LBFGSHessian::evaluate_hessian not implemented");
   }

   void LBFGSHessian::compute_hessian_vector_product(const Model& /*model*/, const Vector<double>& vector, double /*objective_multiplier*/,
         const Vector<double>& /*constraint_multipliers*/, Vector<double>& result) {
      if (this->hessian_recomputation_required) {
         this->compute_hessian_representation();
         this->hessian_recomputation_required = false;
      }

      // for the moment, pretend we have an identity Hessian TODO
      for (size_t variable_index: Range(vector.size())) {
         result[variable_index] = vector[variable_index];
      }
   }

   void LBFGSHessian::update_memory(const Iterate& current_iterate, const Iterate& trial_iterate) {
      std::cout << "Adding vector to L-BFGS memory at slot " << this->current_available_slot << '\n';
      // this->current_available_slot lives in [0, this->memory_size)
      // TODO figure out if we're extending or replacing in memory

      // fill the S matrix
      for (size_t variable_index: Range(this->dimension)) {
         this->S_matrix.entry(variable_index, this->current_available_slot) = trial_iterate.primals[variable_index] -
            current_iterate.primals[variable_index];
      }
      std::cout << "S:\n" << this->S_matrix;

      // fill the Y matrix
      // TODO

      this->used_memory_size = std::min(this->used_memory_size + 1, this->memory_size);
      std::cout << "There are now " << this->used_memory_size << " iterates in memory (capacity " <<
         this->memory_size << ")\n";
   }

   void LBFGSHessian::compute_hessian_representation() {
      // fill the D matrix (diagonal)
      this->D_matrix[this->current_available_slot] = 1.; // TODO dot(s_new, y_new)
      std::cout << "D: "; print_vector(std::cout, this->D_matrix);

      // fill the L matrix (lower triangular with a zero diagonal)
      for (size_t column_index: Range(this->used_memory_size)) {
         for (size_t row_index: Range(column_index+1, this->used_memory_size)) {
            this->L_matrix.entry(row_index, column_index) = 1.; // TODO dot(s_i, y_j)
         }
      }
      std::cout << "L:\n" << this->L_matrix;

      // assemble m x m matrix M = S^T B0 S + L D^{-1} L^T
      // where L D^{-1} L^T = D^{-1} L L^T (because D is diagonal)
      // DenseMatrix<double> Ltilde_matrix(this->current_memory_size, this->current_memory_size);
      // scale
      std::cout << "M:\n" << this->M_matrix;

      // compute the Cholesky factors J of M = J J^T (J overwrites M)
      char uplo = 'L';
      int info = 0;
      int M_dimension = static_cast<int>(this->used_memory_size);
      int M_leading_dimension = static_cast<int>(this->memory_size);
      LAPACK_cholesky_factorization(&uplo, &M_dimension, this->M_matrix.data(), &M_leading_dimension, &info);
      std::cout << "Cholesky info: " << info << '\n';
      std::cout << "J:\n" << this->M_matrix;

      // if we exceed the size of the memory, we start over and replace the older point in memory
      this->current_available_slot = (this->current_available_slot + 1) % this->memory_size;
   }
} // namespace