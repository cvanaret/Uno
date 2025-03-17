// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYQPSOLVER_H
#define UNO_INEQUALITYQPSOLVER_H

#include <vector>

namespace uno {
   // forward declarations
   class Direction;
   //class HessianModel;
   class LagrangeNewtonSubproblem;
   class Statistics;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class InequalityQPSolver {
   public:
      InequalityQPSolver() = default;
      virtual ~InequalityQPSolver() = default;

      virtual void solve_inequality_constrained_QP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) = 0;

      [[nodiscard]] virtual double hessian_quadratic_product(const Vector<double>& primal_direction) const = 0;
   };
} // namespace

#endif // UNO_INEQUALITYQPSOLVER_H