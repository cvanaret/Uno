// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SOLVERWORKSPACE_H
#define UNO_SOLVERWORKSPACE_H

namespace uno {
   // forward declarations
   class Subproblem;
   template <typename ElementType>
   class Vector;

   class SolverWorkspace {
   public:
      SolverWorkspace() = default;
      virtual ~SolverWorkspace() = default;

      [[nodiscard]] virtual double compute_hessian_quadratic_form(const Subproblem& subproblem, const Vector<double>& vector) const = 0;
   };
} // namespace

#endif // UNO_SOLVERWORKSPACE_H