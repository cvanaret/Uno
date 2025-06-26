// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSOLVERFACTORY_H
#define UNO_QPSOLVERFACTORY_H

#include <memory>
#include "QPSolver.hpp"
#include "tools/Literals.hpp"

namespace uno {
   // forward declaration
   class Options;

   class QPSolverFactory {
   public:
      // create a QP solver
      static std::unique_ptr<QPSolver> create([[maybe_unused]] const Options& options);

      // list of available QP solvers
      constexpr static auto available_solvers = make_literal_array(
#ifdef HAS_BQPD
         "BQPD"
#endif
      );
   };
} // namespace

#endif // UNO_QPSOLVERFACTORY_H
