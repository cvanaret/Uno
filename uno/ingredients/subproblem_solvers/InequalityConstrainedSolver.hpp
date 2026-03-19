// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYCONSTRAINEDSOLVER_H
#define UNO_INEQUALITYCONSTRAINEDSOLVER_H

#include "SubproblemSolver.hpp"

namespace uno {
   class InequalityConstrainedSolver: public SubproblemSolver {
   public:
      InequalityConstrainedSolver() = default;
      ~InequalityConstrainedSolver() override = default;
   };
} // namespace

#endif // UNO_INEQUALITYCONSTRAINEDSOLVER_H