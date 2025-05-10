// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include <string>
#include "HessianModelFactory.hpp"
#include "HessianModel.hpp"
#include "ExactHessian.hpp"
#include "IdentityHessian.hpp"
#include "ZeroHessian.hpp"
#include "options/Options.hpp"

namespace uno {
   std::unique_ptr<HessianModel> HessianModelFactory::create(const Options& options) {
      const std::string& hessian_model = options.get_string("hessian_model");
      if (hessian_model == "exact") {
         return std::make_unique<ExactHessian>();
      }
      else if (hessian_model == "identity") {
         return std::make_unique<IdentityHessian>();
      }
      else if (hessian_model == "zero") {
         return std::make_unique<ZeroHessian>();
      }
      throw std::invalid_argument("Hessian model " + hessian_model + " does not exist");
   }
} // namespace