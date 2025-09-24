// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include "InequalityHandlingMethod.hpp"
#include "InequalityHandlingMethodFactory.hpp"
#include "inequality_constrained_methods/InequalityConstrainedMethod.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "interior_point_methods/InteriorPointMethod.hpp"
#include "interior_point_methods/barrier_problems/PrimalExponentialBarrierProblem.hpp"
#include "interior_point_methods/barrier_problems/PrimalDualInteriorPointProblem.hpp"
#include "options/Options.hpp"

namespace uno {
   std::unique_ptr<InequalityHandlingMethod> InequalityHandlingMethodFactory::create(const Options& options) {
      // TODO set unconstrained strategy automatically
      const std::string inequality_handling_method = options.get_string("inequality_handling_method");
      // inequality-constrained methods
      if (inequality_handling_method == "inequality_constrained") {
         return std::make_unique<InequalityConstrainedMethod>(options);
      }
      // interior-point methods
      else if (inequality_handling_method == "interior_point" && options.get_string("barrier_function") == "log") {
         return std::make_unique<InteriorPointMethod<PrimalDualInteriorPointProblem>>(options);
      }
      else if (inequality_handling_method == "interior_point" && options.get_string("barrier_function") == "exp") {
         return std::make_unique<InteriorPointMethod<PrimalExponentialBarrierProblem>>(options);
      }
      throw std::invalid_argument("Inequality handling method " + inequality_handling_method + " is not supported");
   }

   std::vector<std::string> InequalityHandlingMethodFactory::available_strategies() {
      std::vector<std::string> strategies{};
      if constexpr (0 < LPSolverFactory::available_solvers.size() || 0 < QPSolverFactory::available_solvers.size()) {
         strategies.emplace_back("inequality_constrained");
      }
      if (!SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
         strategies.emplace_back("interior_point");
      }
      return strategies;
   }
} // namespace