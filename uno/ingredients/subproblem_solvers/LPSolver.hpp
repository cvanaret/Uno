// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVER_H
#define UNO_LPSOLVER_H

#include "SubproblemStatus.hpp"

namespace uno {
   // forward declarations
   class LagrangeNewtonSubproblem;
   class Multipliers;
   class Statistics;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class LPSolver {
   public:
      LPSolver() = default;
      virtual ~LPSolver() = default;

      [[nodiscard]] virtual SubproblemStatus solve_LP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem,
         const Vector<double>& initial_point, Vector<double>& direction_primals, Multipliers& direction_multipliers, double& subproblem_objective,
         const WarmstartInformation& warmstart_information) = 0;

      [[nodiscard]] double hessian_quadratic_product(const Vector<double>& /*primal_direction*/) const { return 0.; }
   };
} // namespace

#endif // UNO_LPSOLVER_H
