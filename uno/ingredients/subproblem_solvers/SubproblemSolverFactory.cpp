// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "SubproblemSolverFactory.hpp"
#include "BoxLPSolver.hpp"
#include "EQPSolver.hpp"
#include "LPSolverFactory.hpp"
#include "QPSolverFactory.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "tools/Logger.hpp"

namespace uno {
   std::unique_ptr<SubproblemSolver> SubproblemSolverFactory::create(const Subproblem& subproblem, bool uses_trust_region,
         const Options& options) {
      // if no inequality constraint and no trust region, allocate EQP solver
      if (!subproblem.has_inequality_constraints() && !uses_trust_region) {
         DEBUG << "No inequality constraints in the subproblem, allocating an EQP solver\n";
         return std::make_unique<EQPSolver>(options);
      }
      // otherwise, allocate LP/QP solver, depending on the presence of curvature in the subproblem
      if (!subproblem.has_curvature()) {
         if (subproblem.number_constraints == 0) {
            DEBUG << "No curvature and only bound constraints in the subproblem, allocating a box LP solver\n";
            return std::make_unique<BoxLPSolver>();
         }
         else {
            DEBUG << "No curvature in the subproblem, allocating an LP solver\n";
            return LPSolverFactory::create(options);
         }
      }
      else {
         DEBUG << "Curvature in the subproblem, allocating a QP solver\n";
         return QPSolverFactory::create(options);
      }
   }
} // namespace