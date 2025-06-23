// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_SYMMETRICINDEFINITELINEARSOLVER_H

#include <cstddef>

namespace uno {
   // forward declarations
   class Direction;
   class EvaluationSpace;
   class Statistics;
   class Subproblem;
   template <typename ElementType>
   class Vector;
   struct WarmstartInformation;

   template <typename ElementType>
   class SymmetricIndefiniteLinearSolver {
   public:
      SymmetricIndefiniteLinearSolver() = default;
      virtual ~SymmetricIndefiniteLinearSolver() = default;

      virtual void initialize_hessian(const Subproblem& subproblem) = 0;
      virtual void initialize_augmented_system(const Subproblem& subproblem) = 0;

      virtual void solve_indefinite_system(const Vector<double>& matrix_values, const Vector<ElementType>& rhs,
         Vector<ElementType>& result) = 0;
      virtual void solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         const WarmstartInformation& warmstart_information) = 0;

      [[nodiscard]] virtual EvaluationSpace& get_evaluation_space() = 0;
   };
} // namespace

#endif // UNO_SYMMETRICINDEFINITELINEARSOLVER_H