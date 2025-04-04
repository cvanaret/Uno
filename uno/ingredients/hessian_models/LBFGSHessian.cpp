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
   LBFGSHessian::LBFGSHessian(const Model& model): HessianModel("L-BFGS"), model(model) {
   }

   bool LBFGSHessian::has_hessian_operator() const {
      return true;
   }

   bool LBFGSHessian::has_hessian_matrix() const {
      return false;
   }

   bool LBFGSHessian::has_curvature() const {
      return true;
   }

   size_t LBFGSHessian::number_nonzeros() const {
      throw std::runtime_error("LBFGSHessian::number_nonzeros should not be called");
   }

   bool LBFGSHessian::is_positive_definite() const {
      return true;
   }

   void LBFGSHessian::initialize(const Model& /*model*/) {
   }

   void LBFGSHessian::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) const {
   }

   void LBFGSHessian::evaluate_hessian(Statistics& /*statistics*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* /*hessian_values*/) {
      throw std::runtime_error("LBFGSHessian::evaluate_hessian should not be called");
   }

   void LBFGSHessian::compute_hessian_vector_product(const double* /*x*/, const double* /*vector*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* /*result*/) {
      throw std::runtime_error("LBFGSHessian::compute_hessian_vector_product not implemented");
   }
} // namespace