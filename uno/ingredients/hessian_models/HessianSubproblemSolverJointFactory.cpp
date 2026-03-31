// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <utility>
#include <stdexcept>
#include <string>
#include "HessianSubproblemSolverJointFactory.hpp"
#include "HessianModel.hpp"
#include "ExactHessian.hpp"
#include "quasi_newton/LBFGSHessian.hpp"
#include "quasi_newton/LSR1Hessian.hpp"
#include "IdentityHessian.hpp"
#include "ZeroHessian.hpp"
#include "ingredients/subproblem_solvers/SubproblemSolver.hpp"
#include "ingredients/subproblem_solvers/SubproblemSolverFactory.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   std::pair<std::unique_ptr<HessianModel>, std::unique_ptr<SubproblemSolver>> HessianSubproblemSolverJointFactory::create(const Model& model,
         const OptimizationProblem& problem, Iterate& current_iterate, InertiaCorrectionStrategy& inertia_correction_strategy,
         bool uses_trust_region, double objective_multiplier, Options& options) {
      // first look at the problem type. If it is an LP, pick zero Hessian model
      if (model.get_problem_type() == ProblemType::LINEAR) {
         // override user defined option
         options.set_string("hessian_model", "zero", true);
         auto hessian_model = std::make_unique<ZeroHessian>(model.number_variables);
         auto subproblem_solver = SubproblemSolverFactory::create(problem, current_iterate, *hessian_model,
            inertia_correction_strategy, uses_trust_region, options);
         return {std::move(hessian_model), std::move(subproblem_solver)};
      }

      // then look at the option hessian_model
      // from now onwards, the problem is nonlinear (QP or NLP)
      const std::string& hessian_model_type = options.get_string("hessian_model");
      if (hessian_model_type == "exact") {
         if (model.has_hessian_matrix() || model.has_hessian_operator()) {
            auto hessian_model = std::make_unique<ExactHessian>(model);
            auto subproblem_solver = SubproblemSolverFactory::create(problem, current_iterate, *hessian_model,
               inertia_correction_strategy, uses_trust_region, options);
            return {std::move(hessian_model), std::move(subproblem_solver)};
         }
         // no Hessian (matrix or operator) is available: pick an quasi-Newton Hessian (L-BFGS for line search, L-SR1
         // for trust-region methods)
         else if (options.get_string("globalization_mechanism") == "LS") {
            WARNING << "An exact Hessian (matrix or operator) was not provided, setting an L-BFGS Hessian instead\n";
            // override user defined option
            options.set_string("hessian_model", "LBFGS", true);
            auto hessian_model = std::make_unique<LBFGSHessian>(model, objective_multiplier, options);
            auto subproblem_solver = SubproblemSolverFactory::create(problem, current_iterate, *hessian_model,
               inertia_correction_strategy, uses_trust_region, options);
            return {std::move(hessian_model), std::move(subproblem_solver)};
         }
         else {
            WARNING << "An exact Hessian (matrix or operator) was not provided, setting an L-SR1 Hessian instead\n";
            // override user defined option
            options.set_string("hessian_model", "LSR1", true);
            auto hessian_model = std::make_unique<LSR1Hessian>(model, objective_multiplier, options);
            auto subproblem_solver = SubproblemSolverFactory::create(problem, current_iterate, *hessian_model,
               inertia_correction_strategy, uses_trust_region, options);
            return {std::move(hessian_model), std::move(subproblem_solver)};
         }
      }
      else if (hessian_model_type == "LBFGS") {
         auto hessian_model = std::make_unique<LBFGSHessian>(model, objective_multiplier, options);
         auto subproblem_solver = SubproblemSolverFactory::create(problem, current_iterate, *hessian_model,
            inertia_correction_strategy, uses_trust_region, options);
         return {std::move(hessian_model), std::move(subproblem_solver)};
      }
      else if (hessian_model_type == "LSR1") {
         auto hessian_model = std::make_unique<LSR1Hessian>(model, objective_multiplier, options);
         auto subproblem_solver = SubproblemSolverFactory::create(problem, current_iterate, *hessian_model,
            inertia_correction_strategy, uses_trust_region, options);
         return {std::move(hessian_model), std::move(subproblem_solver)};
      }
      else if (hessian_model_type == "identity") {
         auto hessian_model = std::make_unique<IdentityHessian>(model.number_variables);
         auto subproblem_solver = SubproblemSolverFactory::create(problem, current_iterate, *hessian_model,
            inertia_correction_strategy, uses_trust_region, options);
         return {std::move(hessian_model), std::move(subproblem_solver)};
      }
      else if (hessian_model_type == "zero") {
         auto hessian_model = std::make_unique<ZeroHessian>(model.number_variables);
         auto subproblem_solver = SubproblemSolverFactory::create(problem, current_iterate, *hessian_model,
            inertia_correction_strategy, uses_trust_region, options);
         return {std::move(hessian_model), std::move(subproblem_solver)};
      }
      throw std::invalid_argument("Hessian model " + hessian_model_type + " does not exist");
   }
} // namespace