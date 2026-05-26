// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include "InequalityHandlingMethod.hpp"
#include "InequalityHandlingMethodFactory.hpp"
#include "NoInequalityReformulation.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "interior_point_methods/InteriorPointMethod.hpp"
#include "interior_point_methods/barrier_problems/PrimalDualInteriorPointProblem.hpp"
#include "options/Options.hpp"

namespace uno {
   std::unique_ptr<InequalityHandlingMethod> InequalityHandlingMethodFactory::create(const OptimizationProblem& problem,
         bool uses_trust_region, const Options& options) {
      // figure out whether there are inequality constraints altogether
      if (!problem.has_inequality_constraints() && !problem.has_bound_constraints() && !uses_trust_region) {
         // no reformulation
         if (0 < problem.number_constraints) {
            INFO << "The problem has no inequalities, picking a pure SQP method\n";
            return std::make_unique<NoInequalityReformulation>("pure SQP method");
         }
         else {
            INFO << "The problem has no constraints, picking a pure Newton method\n";
            return std::make_unique<NoInequalityReformulation>("pure Newton method");
         }
      }
      // from now on, the problem has inequalities
      const std::string inequality_handling_method = options.get_string("inequality_handling_method");
      // inequality-constrained methods
      if (inequality_handling_method == "inequality_constrained") {
         // no inequality reformulation: let the subproblem solver handle them
         return std::make_unique<NoInequalityReformulation>("inequality-constrained SQP method");
      }
      // interior-point method
      else if (inequality_handling_method == "interior_point") {
         const std::string barrier_function = options.get_string("barrier_function");
         if (barrier_function == "log") {
            return std::make_unique<InteriorPointMethod<PrimalDualInteriorPointProblem>>(options);
         }
         else {
            throw std::invalid_argument("The barrier function " + barrier_function + " is not supported");
         }
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