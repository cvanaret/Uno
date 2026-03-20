// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEMSOLVERFACTORY_H
#define UNO_SUBPROBLEMSOLVERFACTORY_H

#include <memory>
#include "SubproblemSolver.hpp"

namespace uno {
   // forward declarations
   class Options;
   class Subproblem;

   class SubproblemSolverFactory {
   public:
      static std::unique_ptr<SubproblemSolver> create(const Subproblem& subproblem, bool uses_trust_region, const Options& options);
   };
} // namespace

#endif // UNO_SUBPROBLEMSOLVERFACTORY_H