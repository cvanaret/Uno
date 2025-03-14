// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_SYMMETRICINDEFINITELINEARSOLVER_H

#include <cstddef>
#include <tuple>

namespace uno {
   // forward declarations
   class LagrangeNewtonSubproblem;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   template <typename IndexType, typename ElementType>
   class SymmetricIndefiniteLinearSolver {
   public:
      explicit SymmetricIndefiniteLinearSolver(size_t dimension) : dimension(dimension) { };
      virtual ~SymmetricIndefiniteLinearSolver() = default;

      virtual void solve_EQP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, Vector<ElementType>& result,
            WarmstartInformation& warmstart_information) = 0;

   protected:
      const size_t dimension;
   };
} // namespace

#endif // UNO_SYMMETRICINDEFINITELINEARSOLVER_H