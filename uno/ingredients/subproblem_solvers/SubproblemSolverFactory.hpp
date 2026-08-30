// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEMSOLVERFACTORY_H
#define UNO_SUBPROBLEMSOLVERFACTORY_H

#include <memory>
#include "BoxLPSolver.hpp"
#include "EQPSolver.hpp"
#include "IQPSolver.hpp"
#include "LPSolverFactory.hpp"
#include "InverseNewtonSolver.hpp"
#include "QPSolver.hpp"
#include "QPSolverFactory.hpp"
#include "WoodburyEQPSolver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   // forward declarations
   class InverseLBFGSHessian;
   class DirectQuasiNewtonHessian;
   class Options;

   class SubproblemSolverFactory {
   public:
      template <typename HessianType>
      static std::unique_ptr<SubproblemSolver> create(HessianType& hessian_model, const Subproblem& subproblem,
         bool uses_trust_region, const Options& options);
   };

   inline void warn_if_kkt_dump_unsupported(const Options& options) {
      if (!options.get_string("dump_kkt_path").empty()) {
         WARNING << "Uno: 'dump_kkt_path' is set, but KKT system dumping is only supported with an exact Hessian; "
            "no systems will be dumped with the quasi-Newton Hessian\n";
      }
   }

   template <typename HessianType>
   std::unique_ptr<SubproblemSolver> SubproblemSolverFactory::create(HessianType& hessian_model, const Subproblem& subproblem,
         bool uses_trust_region, const Options& options) {
      // if no curvature, allocate LP solver
      if (!subproblem.has_curvature()) {
         if (subproblem.number_constraints == 0) {
            DEBUG << "No curvature and only bound constraints in the subproblem, allocating a box LP solver\n";
            return std::make_unique<BoxLPSolver>();
         }
         else {
            DEBUG << "No curvature in the subproblem, allocating an LP solver\n";
            return std::make_unique<IQPSolver>(LPSolverFactory::create(options));
         }
      }
      // if no inequality constraint and no trust region, allocate EQP solver
      else if (!subproblem.has_inequality_constraints() && !subproblem.has_bound_constraints() && !uses_trust_region) {
         if constexpr (std::is_same_v<HessianType, InverseLBFGSHessian>) { // unconstrained
            DEBUG << "No constraints in the subproblem, allocating a Newton solver with inverse quasi-Newton Hessian\n";
            warn_if_kkt_dump_unsupported(options);
            // the hessian_model we pass has type QuasiNewtonHessian
            return std::make_unique<InverseNewtonSolver>(hessian_model);
         }
         else if constexpr (std::is_base_of_v<DirectQuasiNewtonHessian, HessianType>) { // equality-constrained
            DEBUG << "No inequality constraints in the subproblem, allocating an EQP solver with quasi-Newton Hessian\n";
            warn_if_kkt_dump_unsupported(options);
            // the hessian_model we pass has type QuasiNewtonHessian
            return std::make_unique<WoodburyEQPSolver>(hessian_model, options);
         }
         else {
            DEBUG << "No inequality constraints in the subproblem, allocating an EQP solver\n";
            return std::make_unique<EQPSolver>(options);
         }
      }
      // otherwise, allocate QP solver
      else {
         DEBUG << "Curvature in the subproblem, allocating a QP solver\n";
         // wrap the (Subproblem-agnostic) QP backend in an IQPSolver that builds its QuadraticProgram
         return std::make_unique<IQPSolver>(QPSolverFactory::create(options));
      }
   }
} // namespace

#endif // UNO_SUBPROBLEMSOLVERFACTORY_H
