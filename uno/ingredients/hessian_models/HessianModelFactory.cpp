// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include <string>
#include "HessianModelFactory.hpp"
#include "HessianModel.hpp"
#include "ExactHessian.hpp"
#include "quasi_newton/LBFGSHessian.hpp"
#include "quasi_newton/LSR1Hessian.hpp"
#include "IdentityHessian.hpp"
#include "ZeroHessian.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   std::unique_ptr<HessianModel> HessianModelFactory::create(const Model& model, [[maybe_unused]] double objective_multiplier,
         Options& options) {
      // first look at the problem type. If it is an LP, pick zero Hessian model
      if (model.get_problem_type() == ProblemType::LINEAR) {
         // override user defined option
         options.set_string("hessian_model", "zero", true);
         return std::make_unique<ZeroHessian>(model.number_variables);
      }

      // then look at the option hessian_model
      // from now onwards, the problem is nonlinear (QP or NLP)
      const std::string& hessian_model = options.get_string("hessian_model");
      if (hessian_model == "exact") {
         if (model.has_hessian_matrix() || model.has_hessian_operator()) {
            return std::make_unique<ExactHessian>(model);
         }
         // if no Hessian (matrix or operator) is available, pick a limited-memory quasi-Newton Hessian
         // (L-BFGS for line search, L-SR1 for trust-region method)
         else {
            if (options.get_string("globalization_mechanism") == "LS") {
               WARNING << "An exact Hessian (matrix or operator) was not provided, setting an L-BFGS Hessian instead\n";
               // override user defined option
               options.set_string("hessian_model", "LBFGS", true);
               return std::make_unique<LBFGSHessian>(model, objective_multiplier, options);
            }
            else {
               WARNING << "An exact Hessian (matrix or operator) was not provided, setting an L-SR1 Hessian instead\n";
               // override user defined option
               options.set_string("hessian_model", "LSR1", true);
               return std::make_unique<LSR1Hessian>(model, objective_multiplier, options);
            }
         }
      }
      else if (hessian_model == "LBFGS") {
         return std::make_unique<LBFGSHessian>(model, objective_multiplier, options);
      }
      else if (hessian_model == "LSR1") {
         return std::make_unique<LSR1Hessian>(model, objective_multiplier, options);
      }
      else if (hessian_model == "identity") {
         return std::make_unique<IdentityHessian>(model.number_variables);
      }
      else if (hessian_model == "zero") {
         return std::make_unique<ZeroHessian>(model.number_variables);
      }
      else if (hessian_model == "L-BFGS") {
         return std::make_unique<LBFGSHessian>(model, objective_multiplier, options);
      }
      throw std::invalid_argument("Hessian model " + hessian_model + " does not exist");
   }
} // namespace