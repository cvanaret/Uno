// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LINEARSOLVERFACTORY_H
#define UNO_LINEARSOLVERFACTORY_H

#include <memory>
#include "SymmetricIndefiniteLinearSolver.hpp"

#ifdef HAS_MA57
#include "MA57Solver.hpp"
#endif

class LinearSolverFactory {
public:
   static std::unique_ptr<SymmetricIndefiniteLinearSolver<double>> create(const std::string& linear_solver_name, size_t max_dimension,
         size_t max_number_nonzeros) {
#ifdef HAS_MA57
      if (linear_solver_name == "MA57") {
         return std::make_unique<MA57Solver>(max_dimension, max_number_nonzeros);
      }
#endif
      throw std::invalid_argument("Linear solver name is unknown");
   }
};

#endif // UNO_LINEARSOLVERFACTORY_H
