// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LINEARSOLVERFACTORY_H
#define UNO_LINEARSOLVERFACTORY_H

#include <memory>
#include <vector>
#include "ingredients/subproblem_solvers/EqualityQPSolver.hpp"

namespace uno {
   // forward declarations
   class Options;

   class EqualityQPSolverFactory {
   public:
      static std::unique_ptr<EqualityQPSolver<size_t, double>> create(size_t number_variables, size_t number_constraints,
         size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros, const Options& options);

      // return the list of available solvers
      static std::vector<std::string> available_solvers();
   };
} // namespace

#endif // UNO_LINEARSOLVERFACTORY_H
