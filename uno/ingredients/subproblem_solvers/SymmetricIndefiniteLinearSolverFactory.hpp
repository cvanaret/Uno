// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LINEARSOLVERFACTORY_H
#define UNO_LINEARSOLVERFACTORY_H

#include <memory>
#include <string>
#include <vector>

namespace uno {
   // forward declaration
   template <class ElementType>
   class DirectSymmetricIndefiniteLinearSolver;

   class SymmetricIndefiniteLinearSolverFactory {
   public:
      static std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> create(const std::string& linear_solver);

      // return the list of available solvers
      static std::vector<std::string> available_solvers();
   };
} // namespace

#endif // UNO_LINEARSOLVERFACTORY_H