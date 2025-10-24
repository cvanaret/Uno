// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include <string>
#include "HessianModelFactory.hpp"
#include "HessianModel.hpp"
#include "ExactHessian.hpp"
#include "IdentityHessian.hpp"
#include "ZeroHessian.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   std::unique_ptr<HessianModel> HessianModelFactory::create(const Model& model, const Options& options) {
      const std::string& hessian_model = options.get_string("hessian_model");
      if (hessian_model == "exact") {
         // if no Hessian (matrix or operator) is available, pick a zero Hessian
         if (!model.has_hessian_matrix() && !model.has_hessian_operator()) {
            WARNING << "An exact Hessian (matrix or operator) was not provided, setting a zero Hessian instead\n";
            return std::make_unique<ZeroHessian>(model.number_variables);
         }
         else {
            return std::make_unique<ExactHessian>(model);
         }
      }
      else if (hessian_model == "identity") {
         return std::make_unique<IdentityHessian>(model.number_variables);
      }
      else if (hessian_model == "zero") {
         return std::make_unique<ZeroHessian>(model.number_variables);
      }
      throw std::invalid_argument("Hessian model " + hessian_model + " does not exist");
   }
} // namespace