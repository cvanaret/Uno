// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "HessianModelFactory.hpp"
#include "HessianModel.hpp"
#include "ConvexifiedHessian.hpp"
#include "ExactHessian.hpp"
#include "LBFGSHessian.hpp"
#include "ZeroHessian.hpp"

namespace uno {
   std::unique_ptr<HessianModel> HessianModelFactory::create(const std::string& hessian_model, size_t dimension, size_t number_hessian_nonzeros,
         bool convexify, const Options& options) {
      if (hessian_model == "exact") {
         if (convexify) {
            return std::make_unique<ConvexifiedHessian>(dimension, number_hessian_nonzeros, options);
         }
         else {
            return std::make_unique<ExactHessian>();
         }
      }
#ifdef HAS_LAPACK
      else if (hessian_model == "L-BFGS") {
         return std::make_unique<LBFGSHessian>(dimension);
      }
#endif
      else if (hessian_model == "zero") {
         return std::make_unique<ZeroHessian>();
      }
      throw std::invalid_argument("Hessian model " + hessian_model + " does not exist");
   }
} // namespace