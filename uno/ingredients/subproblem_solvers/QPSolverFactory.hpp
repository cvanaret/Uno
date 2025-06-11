// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSOLVERFACTORY_H
#define UNO_QPSOLVERFACTORY_H

#include <initializer_list>
#include <memory>
#include "QPSolver.hpp"

namespace uno {
   // forward declaration
   class Options;

   class QPSolverFactory {
   public:
      // create a QP solver
      static std::unique_ptr<QPSolver> create([[maybe_unused]] const Options& options);

      // list of available QP solvers
      constexpr static std::initializer_list<const char*> available_solvers{
#ifdef HAS_BQPD
         "BQPD",
#endif
#ifdef HAS_HIGHS
         "HiGHS",
#endif
      };
   };
} // namespace

#endif // UNO_QPSOLVERFACTORY_H