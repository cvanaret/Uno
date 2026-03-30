// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEMSOLVERFACTORY_H
#define UNO_SUBPROBLEMSOLVERFACTORY_H

#include <memory>
#include "BoxLPSolver.hpp"
#include "EQPSolver.hpp"
#include "LPSolver.hpp"
#include "LPSolverFactory.hpp"
#include "QPSolver.hpp"
#include "QPSolverFactory.hpp"
#include "WoodburyEQPSolver.hpp"
#include "dogleg/DoglegMethod.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   // forward declarations
   class InertiaCorrectionStrategy;
   class LBFGSHessian;
   class Options;
   class Subproblem;

   class SubproblemSolverFactory {
   public:
      template <typename HessianType>
      static std::unique_ptr<SubproblemSolver> create(const OptimizationProblem& problem, Iterate& current_iterate,
         HessianType& hessian_model, InertiaCorrectionStrategy& inertia_correction_strategy,
         bool uses_trust_region, const Options& options);
   };

   template <typename HessianType>
   inline std::unique_ptr<SubproblemSolver> SubproblemSolverFactory::create(const OptimizationProblem& problem,
         Iterate& current_iterate, HessianType& hessian_model, InertiaCorrectionStrategy& inertia_correction_strategy,
         bool uses_trust_region, const Options& options) {
      const Subproblem subproblem(problem, current_iterate, hessian_model, inertia_correction_strategy);
      // if no inequality constraint and no trust region, allocate EQP solver
      // temporary fix: this is set only in interior-point methods
      /*if (!subproblem.has_inequality_constraints() && uses_trust_region && subproblem.is_hessian_positive_definite()) {
         // use the dogleg method
         DEBUG << "Trust-region and no inequality constraints in the subproblem, allocating a dogleg solver\n";
         auto subproblem_solver = std::make_unique<DoglegMethod>(options);
         subproblem_solver->initialize_memory(subproblem);
         return subproblem_solver;
      }
      else*/
      if (!subproblem.has_inequality_constraints() && !uses_trust_region && options.get_string("inequality_handling_method") == "interior_point") {
         // no trust region
         if constexpr (std::is_same_v<HessianType, LBFGSHessian>) {
            DEBUG << "No inequality constraints in the subproblem, allocating an EQP solver with L-BFGS Hessian\n";
            // the hessian_model we pass has type LBFGSHessian
            auto subproblem_solver = std::make_unique<WoodburyEQPSolver>(hessian_model, options);
            subproblem_solver->initialize_memory(subproblem);
            return subproblem_solver;
         }
         else {
            DEBUG << "No inequality constraints in the subproblem, allocating an EQP solver\n";
            auto subproblem_solver = std::make_unique<EQPSolver>(options);
            subproblem_solver->initialize_memory(subproblem);
            return subproblem_solver;
         }
      }
      // otherwise, allocate LP/QP solver, depending on the presence of curvature in the subproblem
      else if (!subproblem.has_curvature()) {
         if (subproblem.number_constraints == 0) {
            DEBUG << "No curvature and only bound constraints in the subproblem, allocating a box LP solver\n";
            auto subproblem_solver = std::make_unique<BoxLPSolver>();
            subproblem_solver->initialize_memory(subproblem);
            return subproblem_solver;
         }
         else {
            DEBUG << "No curvature in the subproblem, allocating an LP solver\n";
            auto subproblem_solver = LPSolverFactory::create(options);
            subproblem_solver->initialize_memory(subproblem);
            return subproblem_solver;
         }
      }
      else {
         DEBUG << "Curvature in the subproblem, allocating a QP solver\n";
         auto subproblem_solver = QPSolverFactory::create(options);
         subproblem_solver->initialize_memory(subproblem);
         return subproblem_solver;
      }
   }
} // namespace

#endif // UNO_SUBPROBLEMSOLVERFACTORY_H