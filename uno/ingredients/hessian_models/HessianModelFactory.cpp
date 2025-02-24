// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "HessianModelFactory.hpp"
#include "HessianModel.hpp"
#include "ConvexifiedHessian.hpp"
#include "ExactHessian.hpp"
#include "IdentityHessian.hpp"
#include "ZeroHessian.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"

namespace uno {
   std::unique_ptr<HessianModel> HessianModelFactory::create(const std::string& hessian_model, size_t dimension, size_t maximum_number_nonzeros,
         bool regularize, const Options& options) {
      if (hessian_model == "exact") {
         if (regularize) {
            return std::make_unique<ConvexifiedHessian>(dimension, maximum_number_nonzeros + dimension, options);
         }
         else {
            return std::make_unique<ExactHessian>();
         }
      }
      else if (hessian_model == "zero") {
         return std::make_unique<ZeroHessian>();
      }
      else if (hessian_model == "identity") {
         return std::make_unique<IdentityHessian>();
      }
      throw std::invalid_argument("Hessian model " + hessian_model + " does not exist");
   }
} // namespace