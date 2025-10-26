// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "fortran_interface.h"
#include "LBFGSHessian.hpp"

#define LAPACK_cholesky_factorization FC_GLOBAL_(dpotrf, DPOTRF)

extern "C" {
   void LAPACK_cholesky_factorization(char* uplo, int* n, double* a, int* lda, int* info);
}

namespace uno {
   LBFGSHessian::LBFGSHessian(): HessianModel() {
   }

   bool LBFGSHessian::has_implicit_representation() const {
      return true;
   }

   bool LBFGSHessian::has_explicit_representation() const {
      return false;
   }

   void LBFGSHessian::initialize(const Model& /*model*/) {
      // do nothing
   }

   size_t LBFGSHessian::number_nonzeros(const Model& /*model*/) const {
      throw std::runtime_error("LBFGSHessian::number_nonzeros should not be called");
   }

   bool LBFGSHessian::is_positive_definite() const {
      return true;
   }

   void LBFGSHessian::evaluate_hessian(Statistics& /*statistics*/, const Model& /*model*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, SymmetricMatrix<size_t, double>& /*hessian*/) {
      throw std::runtime_error("LBFGSHessian::evaluate_hessian should not be called");
   }

   void LBFGSHessian::compute_hessian_vector_product(const Model& /*model*/, const double* /*vector*/, double /*objective_multiplier*/,
         const Vector<double>& /*constraint_multipliers*/, double* /*result*/) {
      throw std::runtime_error("LBFGSHessian::compute_hessian_vector_product not implemented");
   }

   std::string LBFGSHessian::get_name() const {
      return "L-BFGS";
   }
} // namespace