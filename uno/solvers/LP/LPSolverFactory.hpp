// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVERFACTORY_H
#define UNO_LPSOLVERFACTORY_H

#include <memory>
#include "LPSolver.hpp"

#ifdef HAS_BQPD
#include "solvers/QP/BQPDSolver.hpp"
#endif

class LPSolverFactory {
public:
   static std::unique_ptr<LPSolver> create(const std::string& LP_solver_name, size_t number_variables, size_t number_constraints,
         size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros, const Options& options) {
#ifdef HAS_BQPD
      if (LP_solver_name == "BQPD") {
         return std::make_unique<BQPDSolver>(number_variables, number_constraints, number_objective_gradient_nonzeros,
               number_jacobian_nonzeros, 0, BQPDProblemType::LP, options);
      }
#endif
      throw std::invalid_argument("LP solver not found");
   }
};

#endif // UNO_LPSOLVERFACTORY_H
