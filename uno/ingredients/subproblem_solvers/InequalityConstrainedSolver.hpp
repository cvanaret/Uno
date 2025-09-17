// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYCONSTRAINEDSOLVER_H
#define UNO_INEQUALITYCONSTRAINEDSOLVER_H

namespace uno {
   // forward declarations
   class Direction;
   class EvaluationSpace;
   class Statistics;
   class Subproblem;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class InequalityConstrainedSolver {
   public:
      InequalityConstrainedSolver() = default;
      virtual ~InequalityConstrainedSolver() = default;

      virtual void initialize_memory(const Subproblem& subproblem) = 0;

      virtual void solve(Statistics& statistics, Subproblem& subproblem, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) = 0;

      [[nodiscard]] virtual EvaluationSpace& get_evaluation_space() = 0;
   };
} // namespace

#endif // UNO_INEQUALITYCONSTRAINEDSOLVER_H