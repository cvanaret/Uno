// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSOLVERFACTORY_H
#define UNO_QPSOLVERFACTORY_H

#include <memory>
#include <vector>

namespace uno {
   // forward declarations
   class Options;
   class QPSolver;

   class QPSolverFactory {
   public:
      // create a QP solver
      static std::unique_ptr<QPSolver> create([[maybe_unused]] const Options& options);

      // return the list of available QP solvers
      static std::vector<std::string> available_solvers();
   };
} // namespace

#endif // UNO_QPSOLVERFACTORY_H
