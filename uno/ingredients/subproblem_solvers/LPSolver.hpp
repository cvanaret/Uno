// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVER_H
#define UNO_LPSOLVER_H

namespace uno {
   // forward declarations
   class Direction;
   class HessianModel;
   class Iterate;
   class Multipliers;
   class OptimizationProblem;
   template <typename ElementType>
   class RegularizationStrategy;
   class Statistics;
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

      virtual void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, const Vector<double>& initial_point, Direction& direction,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius,
         const WarmstartInformation& warmstart_information) = 0;

      [[nodiscard]] virtual double hessian_quadratic_product(const Vector<double>& vector) const = 0;
   };
} // namespace

#endif // UNO_LPSOLVER_H
