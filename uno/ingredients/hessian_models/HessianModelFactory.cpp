// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "HessianModelFactory.hpp"
#include "ExactHessian.hpp"
#include "IdentityHessian.hpp"
#include "ZeroHessian.hpp"

namespace uno {
   // forward declaration
   class HessianModel;

   std::unique_ptr<HessianModel> HessianModelFactory::create(const std::string& hessian_model) {
      if (hessian_model == "exact") {
         return std::make_unique<ExactHessian>();
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