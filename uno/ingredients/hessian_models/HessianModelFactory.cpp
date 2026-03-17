// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include <string>
#include "HessianModelFactory.hpp"
#include "HessianModel.hpp"
#include "ExactHessian.hpp"
#include "quasi_newton/LBFGSHessian.hpp"
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
         // no Hessian (matrix or operator) is available
         // SQP methods: pick an L-BFGS Hessian
         else if (options.get_string("inequality_handling_method") == "inequality_constrained") {
            WARNING << "An exact Hessian (matrix or operator) was not provided, setting an L-BFGS Hessian instead\n";
            // override user defined option
            options.set_string("hessian_model", "LBFGS", true);
            return std::make_unique<LBFGSHessian>(model, objective_multiplier, options);
         }
         else { // interior-point methods: pick an identity Hessian for now
            WARNING << "An exact Hessian (matrix or operator) was not provided, setting an identity Hessian instead\n";
            // override user defined option
            options.set_string("hessian_model", "identity", true);
            return std::make_unique<IdentityHessian>(model.number_variables);
         }
      }
      else if (hessian_model == "LBFGS") {
         // SQP methods: pick an L-BFGS Hessian
         if (options.get_string("inequality_handling_method") == "inequality_constrained") {
            return std::make_unique<LBFGSHessian>(model, objective_multiplier, options);
         }
         else { // interior-point methods: pick an identity Hessian for now
            WARNING << "An exact Hessian (matrix or operator) was not provided, setting an identity Hessian instead\n";
            // override user defined option
            options.set_string("hessian_model", "identity", true);
            return std::make_unique<IdentityHessian>(model.number_variables);
         }
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