// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LINEARSOLVERFACTORY_H
#define UNO_LINEARSOLVERFACTORY_H

#include <memory>
#include <vector>

namespace uno {
   // forward declarations
   template <class IndexType, class ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class Options;

   class SymmetricIndefiniteLinearSolverFactory {
   public:
      static std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<size_t, double>> create([[maybe_unused]] size_t dimension,
            [[maybe_unused]] size_t number_nonzeros, const Options& options);

      // return the list of available solvers
      static std::vector<std::string> available_solvers();
   };
} // namespace

#endif // UNO_LINEARSOLVERFACTORY_H
