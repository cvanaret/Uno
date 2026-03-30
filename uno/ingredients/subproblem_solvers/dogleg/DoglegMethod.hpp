// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DOGLEGMETHOD_H
#define UNO_DOGLEGMETHOD_H

#include "../SubproblemSolver.hpp"
#include "DoglegWorkspace.hpp"

namespace uno {
   // forward declaration
   class Options;

   class DoglegMethod: public SubproblemSolver {
   public:
      explicit DoglegMethod(const Options& options);
      ~DoglegMethod() override = default;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& initial_point, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   protected:
      DoglegWorkspace workspace;
   };
} // namespace

#endif // UNO_DOGLEGMETHOD_H