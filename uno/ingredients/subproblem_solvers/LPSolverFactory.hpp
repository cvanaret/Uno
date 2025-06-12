// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVERFACTORY_H
#define UNO_LPSOLVERFACTORY_H

#include <array>
#include <memory>

namespace uno {
   // forward declarations
   class LPSolver;
   class Options;

   class LPSolverFactory {
   public:
      static std::unique_ptr<LPSolver> create([[maybe_unused]] const Options& options);

      // list of available LP solvers
      constexpr static std::array available_solvers{
#ifdef HAS_BQPD
         "BQPD",
#endif
#ifdef HAS_HIGHS
         "HiGHS"
#endif
      };
   };
} // namespace

#endif // UNO_LPSOLVERFACTORY_H