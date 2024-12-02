// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVERFACTORY_H
#define UNO_LPSOLVERFACTORY_H

#include <memory>
#include <vector>

namespace uno {
   // forward declarations
   class LPSolver;
   class Options;

   class LPSolverFactory {
   public:
      static std::unique_ptr<LPSolver> create([[maybe_unused]] size_t number_variables, [[maybe_unused]] size_t number_constraints,
            [[maybe_unused]] size_t number_objective_gradient_nonzeros, [[maybe_unused]] size_t number_jacobian_nonzeros,
            [[maybe_unused]] const Options& options);

      static std::vector<std::string> available_solvers();
   };
} // namespace

#endif // UNO_LPSOLVERFACTORY_H
