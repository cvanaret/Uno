// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_SYMMETRICINDEFINITELINEARSOLVER_H

#include <cstddef>
#include <tuple>
#include "SubproblemStatus.hpp"

namespace uno {
   // forward declarations
   class LagrangeNewtonSubproblem;
   class Multipliers;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   template <typename IndexType, typename ElementType>
   class EqualityQPSolver {
   public:
      EqualityQPSolver() = default;
      virtual ~EqualityQPSolver() = default;

      virtual SubproblemStatus solve_equality_constrained_QP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem,
         const Vector<double>& initial_point, Vector<double>& direction_primals, Multipliers& direction_multipliers, double& subproblem_objective,
         WarmstartInformation& warmstart_information) = 0;
   };
} // namespace

#endif // UNO_SYMMETRICINDEFINITELINEARSOLVER_H