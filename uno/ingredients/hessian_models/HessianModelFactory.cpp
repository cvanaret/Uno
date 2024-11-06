// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HessianModelFactory.hpp"
#include "HessianModel.hpp"
#include "ConvexifiedHessian.hpp"
#include "ExactHessian.hpp"
#include "ZeroHessian.hpp"
#include "solvers/DirectSymmetricIndefiniteLinearSolver.hpp"

namespace uno {
   std::unique_ptr<HessianModel> HessianModelFactory::create(const std::string& hessian_model, size_t dimension, size_t maximum_number_nonzeros,
         bool convexify, const Options& options) {
      if (hessian_model == "exact") {
         if (convexify) {
            return std::make_unique<ConvexifiedHessian>(dimension, maximum_number_nonzeros + dimension, options);
         }
         else {
            return std::make_unique<ExactHessian>(dimension, maximum_number_nonzeros, options);
         }
      }
      else if (hessian_model == "zero") {
         return std::make_unique<ZeroHessian>(dimension, options);
      }
      throw std::invalid_argument("Hessian model " + hessian_model + " does not exist");
   }
} // namespace