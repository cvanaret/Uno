// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "LBFGSHessian.hpp"

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
         memory_size(memory_size),
         S_matrix(dimension, memory_size),
         Y_matrix(dimension, memory_size),
         M_matrix(memory_size, memory_size) {
   }

   void LBFGSHessian::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) const {
      // do nothing
   }

   void LBFGSHessian::notify_accepted_iterate(const Iterate& /*iterate*/) {

   }

   void LBFGSHessian::evaluate_hessian(Statistics& /*statistics*/, const OptimizationProblem& /*problem*/,
         const Vector<double>& /*primal_variables*/, const Vector<double>& /*constraint_multipliers*/,
         SymmetricMatrix<size_t, double>& /*hessian*/) {
      throw std::runtime_error("LBFGSHessian::evaluate_hessian not implemented");
   }

   void LBFGSHessian::compute_hessian_vector_product(const OptimizationProblem& /*problem*/, const Vector<double>& /*vector*/,
         const Vector<double>& /*constraint_multipliers*/, Vector<double>& /*result*/) {
      throw std::runtime_error("LBFGSHessian::compute_hessian_vector_product not implemented");
   }
} // namespace