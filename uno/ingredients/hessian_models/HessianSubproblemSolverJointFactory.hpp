// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HESSIANSUBPROBLEMSOLVERJOINTFACTORY_H
#define UNO_HESSIANSUBPROBLEMSOLVERJOINTFACTORY_H

#include <array>
#include <memory>
#include <utility>

namespace uno {
   // forward declarations
   class HessianModel;
   class InertiaCorrectionStrategy;
   class Iterate;
   class Model;
   class OptimizationProblem;
   class Options;
   class SubproblemSolver;

   class HessianSubproblemSolverJointFactory {
   public:
      // joint factory of Hessian models and subproblem solvers
      static std::pair<std::unique_ptr<HessianModel>, std::unique_ptr<SubproblemSolver>> create(const Model& model,
         const OptimizationProblem& problem, Iterate& current_iterate, InertiaCorrectionStrategy& inertia_correction_strategy,
         bool uses_trust_region, double objective_multiplier, Options& options);

      constexpr static std::array available_strategies{"exact", "LFBGS", "LSR1", "identity", "zero"};
   };
} // namespace

#endif // UNO_HESSIANSUBPROBLEMSOLVERJOINTFACTORY_H,