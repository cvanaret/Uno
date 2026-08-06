// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HESSIANSUBPROBLEMSOLVERJOINTFACTORY_H
#define UNO_HESSIANSUBPROBLEMSOLVERJOINTFACTORY_H

#include <array>
#include <memory>
#include <tuple>

namespace uno {
   // forward declarations
   class HessianModel;
   class InertiaCorrectionStrategy;
   class OptimizationProblem;
   class Options;
   class SubproblemSolver;

   using Ingredients = std::tuple<std::unique_ptr<InertiaCorrectionStrategy>,
                                  std::unique_ptr<HessianModel>,
                                  std::unique_ptr<SubproblemSolver>
                                 >;

   class HessianSubproblemSolverJointFactory {
   public:
      // joint factory of inertia correction strategy, Hessian models, and subproblem solver
      static Ingredients create(const OptimizationProblem& problem, bool uses_trust_region, double objective_multiplier,
         Options& options);

      constexpr static std::array available_strategies{"exact", "LFBGS", "LSR1", "identity", "zero"};
   };
} // namespace

#endif // UNO_HESSIANSUBPROBLEMSOLVERJOINTFACTORY_H,