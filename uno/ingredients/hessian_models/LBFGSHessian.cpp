// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "LBFGSHessian.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Range.hpp"

#ifdef WITH_LAPACK
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
         memory_capacity(memory_size),
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
      std::cout << "Adding vector to L-BFGS memory at slot " << this->current_available_slot << '\n';
      // this->current_available_slot lives in [0, this->memory_capacity)

      // TODO figure out if we're extending or replacing in memory

      // fill the S matrix
      for (size_t variable_index: Range(this->dimension)) {
         this->S_matrix.entry(variable_index, this->current_available_slot) = trial_iterate.primals[variable_index] -
            current_iterate.primals[variable_index];
      }
      std::cout << "S:\n" << this->S_matrix;

      // fill the Y matrix
      // TODO

      // fill the D matrix (diagonal)
      this->D_matrix[this->current_available_slot] = 1.; // TODO dot(s_new, y_new)

      // fill the L matrix (lower triangular with a zero diagonal)
      for (size_t column_index: Range(this->current_memory_size)) {
         for (size_t row_index: Range(column_index+1, this->current_memory_size)) {
            this->L_matrix.entry(row_index, column_index) = 1.; // TODO dot(s_i, y_j)
         }
      }
      std::cout << "L:\n" << this->L_matrix;

      // if we exceed the size of the memory, we start over and replace the older point in memory
      this->current_available_slot = (this->current_available_slot + 1) % this->memory_capacity;
      this->current_memory_size = std::min(current_memory_size + 1, this->memory_capacity);
      std::cout << "There are now " << this->current_memory_size << " iterates in memory (capacity " <<
         this->memory_capacity << ")\n";
   }

   void LBFGSHessian::evaluate_hessian(Statistics& /*statistics*/, const Model& /*model*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, SymmetricMatrix<size_t, double>& /*hessian*/) {
      throw std::runtime_error("LBFGSHessian::evaluate_hessian not implemented");
   }

   void LBFGSHessian::compute_hessian_vector_product(const Model& /*model*/, const Vector<double>& vector, double /*objective_multiplier*/,
         const Vector<double>& /*constraint_multipliers*/, Vector<double>& result) {
      //throw std::runtime_error("LBFGSHessian::compute_hessian_vector_product not implemented");
      for (size_t variable_index: Range(vector.size())) {
         result[variable_index] = vector[variable_index];
      }
   }
} // namespace