// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYQPSOLVER_H
#define UNO_INEQUALITYQPSOLVER_H

#include "SubproblemStatus.hpp"

namespace uno {
   // forward declarations
   class Direction;
   //class HessianModel;
   class LagrangeNewtonSubproblem;
   class Multipliers;
   class Statistics;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class InequalityQPSolver {
   public:
      InequalityQPSolver() = default;
      virtual ~InequalityQPSolver() = default;

      virtual SubproblemStatus solve_inequality_constrained_QP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem,
         const Vector<double>& initial_point, Vector<double>& direction_primals, Multipliers& direction_multipliers, double& subproblem_objective,
         const WarmstartInformation& warmstart_information) = 0;

      [[nodiscard]] virtual double hessian_quadratic_product(const Vector<double>& primal_direction) const = 0;
   };
} // namespace

#endif // UNO_INEQUALITYQPSOLVER_H