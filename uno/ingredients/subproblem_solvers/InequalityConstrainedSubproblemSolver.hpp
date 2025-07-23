// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYCONSTRAINEDSUBPROBLEMSOLVER_H
#define UNO_INEQUALITYCONSTRAINEDSUBPROBLEMSOLVER_H

namespace uno {
   // forward declarations
   class Direction;
   class HessianModel;
   class OptimizationProblem;
   template <typename ElementType>
   class RegularizationStrategy;
   class Statistics;
   class Subproblem;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class InequalityConstrainedSubproblemSolver {
   public:
      InequalityConstrainedSubproblemSolver() = default;
      virtual ~InequalityConstrainedSubproblemSolver() = default;

      virtual void initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model,
         const RegularizationStrategy<double>& regularization_strategy) = 0;

      virtual void solve_inequality_constrained_subproblem(Statistics& statistics, Subproblem& subproblem,
         const Vector<double>& initial_point, Direction& direction, const WarmstartInformation& warmstart_information) = 0;

      [[nodiscard]] virtual double hessian_quadratic_product(const Vector<double>& vector) const = 0;
   };
} // namespace

#endif // UNO_INEQUALITYCONSTRAINEDSUBPROBLEMSOLVER_H