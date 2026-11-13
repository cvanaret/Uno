// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSOLVER_H
#define UNO_QPSOLVER_H

#include "LPSolver.hpp"

namespace uno {
   // Numerical backend that additionally solves QuadraticPrograms with a quadratic term (a Hessian).
   // The inheritance encodes "every QP solver is also an LP solver, but not vice versa": a QPSolver*
   // is usable wherever an LPSolver* is expected (so a QP solver may be returned by LPSolverFactory),
   // while an LPSolver is not necessarily a QPSolver. The Hessian travels inside the QuadraticProgram,
   // so this refinement adds no method; it is a capability marker.
   class QPSolver : public LPSolver {
   public:
      QPSolver() = default;
      ~QPSolver() override;
   };
} // namespace

#endif // UNO_QPSOLVER_H
