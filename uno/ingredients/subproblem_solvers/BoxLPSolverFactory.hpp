// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BOXLPSOLVERFACTORY_H
#define UNO_BOXLPSOLVERFACTORY_H

#include <initializer_list>
#include <memory>

namespace uno {
   // forward declaration
   class InequalityConstrainedSolver;

   class BoxLPSolverFactory {
   public:
      static std::unique_ptr<InequalityConstrainedSolver> create();

      // list of available LP solvers
      constexpr static std::initializer_list<const char*> available_solvers{"BoxLPSolver"};
   };
} // namespace

#endif // UNO_BOXLPSOLVERFACTORY_H