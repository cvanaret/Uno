// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYQPSOLVERFACTORY_H
#define UNO_INEQUALITYQPSOLVERFACTORY_H

#include <memory>
#include <vector>

namespace uno {
   // forward declarations
   class Options;
   class InequalityQPSolver;

   class InequalityQPSolverFactory {
   public:
      static std::unique_ptr<InequalityQPSolver> create([[maybe_unused]] size_t number_variables, [[maybe_unused]] size_t number_constraints,
            [[maybe_unused]] size_t number_objective_gradient_nonzeros, [[maybe_unused]] size_t number_jacobian_nonzeros,
            [[maybe_unused]] size_t number_hessian_nonzeros, [[maybe_unused]] const Options& options);

      // return the list of available inequality QP solvers
      static std::vector<std::string> available_solvers();
   };
} // namespace

#endif // UNO_INEQUALITYQPSOLVERFACTORY_H
