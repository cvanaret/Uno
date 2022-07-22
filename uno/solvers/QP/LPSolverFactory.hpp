// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVERFACTORY_H
#define UNO_LPSOLVERFACTORY_H

#include <memory>
#include "LPSolver.hpp"

#ifdef HAS_BQPD
#include "BQPDSolver.hpp"
#endif

class LPSolverFactory {
public:
   static std::unique_ptr<BQPDSolver> create(size_t number_variables, size_t number_constraints, const std::string& LP_solver_name,
         const Options& options) {
#ifdef HAS_BQPD
      if (LP_solver_name == "BQPD") {
         return std::make_unique<BQPDSolver>(number_variables, number_constraints, 0, false, options);
      }
      throw std::invalid_argument("LP solver not found");
#endif
   }
};

#endif // UNO_LPSOLVERFACTORY_H
