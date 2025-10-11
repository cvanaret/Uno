// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include <string>
#include "HessianModelFactory.hpp"
#include "HessianModel.hpp"
#include "ExactHessian.hpp"
#include "IdentityHessian.hpp"
#include "ZeroHessian.hpp"
#ifdef HAS_LAPACK
#include "quasi_newton/LBFGSHessian.hpp"
#endif
#include "options/Options.hpp"

namespace uno {
   std::unique_ptr<HessianModel> HessianModelFactory::create([[maybe_unused]] double objective_multiplier, const Options& options) {
      const std::string& hessian_model = options.get_string("hessian_model");
      if (hessian_model == "exact") {
         return std::make_unique<ExactHessian>();
      }
      else if (hessian_model == "identity") {
         return std::make_unique<IdentityHessian>();
      }
#ifdef HAS_LAPACK
      else if (hessian_model == "L-BFGS") {
         return std::make_unique<LBFGSHessian>(objective_multiplier, options);
      }
#endif
      else if (hessian_model == "zero") {
         return std::make_unique<ZeroHessian>();
      }
      throw std::invalid_argument("Hessian model " + hessian_model + " does not exist");
   }
} // namespace