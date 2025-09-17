// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "BoxLPSolverFactory.hpp"
#include "BoxLPSolver.hpp"
#include "InequalityConstrainedSolver.hpp"

namespace uno {
   std::unique_ptr<InequalityConstrainedSolver> BoxLPSolverFactory::create() {
      return std::make_unique<BoxLPSolver>();
   }
} // namespace