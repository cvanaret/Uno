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
   class Options;

   class SymmetricIndefiniteLinearSolverFactory {
   public:
      // libhsl_path: runtime HSL library name for MA27/MA57 (empty = platform default libhsl.<ext>)
      // allow_scaling: attach the linear system scaling strategy (option linear_system_scaling). Solvers that only compute
      // inertias do not need it (a positive diagonal scaling preserves the inertia)
      static std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> create(const Options& options, bool allow_scaling = true);

      // return the list of available solvers
      static std::vector<std::string> available_solvers();
   };
} // namespace

#endif // UNO_LINEARSOLVERFACTORY_H
