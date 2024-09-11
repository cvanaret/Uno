// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSOLVERFACTORY_H
#define UNO_QPSOLVERFACTORY_H

#include <memory>
#include "QPSolver.hpp"

#ifdef HAS_BQPD
#include "BQPDSolver.hpp"
#endif

namespace uno {
   class QPSolverFactory {
   public:
      // create a QP solver
      static std::unique_ptr<QPSolver> create([[maybe_unused]] const std::string& QP_solver_name, [[maybe_unused]] size_t number_variables,
            [[maybe_unused]] size_t number_constraints, [[maybe_unused]] size_t number_objective_gradient_nonzeros,
            [[maybe_unused]] size_t number_jacobian_nonzeros, [[maybe_unused]] size_t number_hessian_nonzeros, [[maybe_unused]] const Options& options) {
   #ifdef HAS_BQPD
         if (QP_solver_name == "BQPD") {
            return std::make_unique<BQPDSolver>(number_variables, number_constraints, number_objective_gradient_nonzeros, number_jacobian_nonzeros,
               number_hessian_nonzeros, BQPDProblemType::QP, options);
         }
   #endif
         throw std::invalid_argument("QP solver name is unknown");
      }

      // return the list of available QP solvers
      static std::vector<std::string> available_solvers() {
         std::vector<std::string> solvers{};
         #ifdef HAS_BQPD
         solvers.emplace_back("BQPD");
         #endif
         return solvers;
      }
   };
} // namespace

#endif // UNO_QPSOLVERFACTORY_H
