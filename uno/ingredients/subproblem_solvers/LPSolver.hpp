// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVER_H
#define UNO_LPSOLVER_H

namespace uno {
   // forward declarations
   class Direction;
   class HessianModel;
   class Iterate;
   class OptimizationProblem;
   template <typename ElementType>
   class RegularizationStrategy;
   template <typename ElementType>
   class Vector;
   struct WarmstartInformation;

   /*! \class LPSolver
    * \brief LP solver
    *
    */
   class LPSolver {
   public:
      LPSolver() = default;
      virtual ~LPSolver() = default;

      virtual void initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model,
         RegularizationStrategy<double>& regularization_strategy) = 0;

      virtual void solve_LP(const OptimizationProblem& problem, Iterate& current_iterate, const Vector<double>& initial_point, Direction& direction,
            double trust_region_radius, const WarmstartInformation& warmstart_information) = 0;
   };
} // namespace

#endif // UNO_LPSOLVER_H
