// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DOGLEGMETHOD_H
#define UNO_DOGLEGMETHOD_H

#include <memory>
#include <string>
#include "../SubproblemSolver.hpp"
#include "DoglegWorkspace.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"

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
      const std::string& linear_solver_name;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
      DoglegWorkspace workspace{};
   };
} // namespace

#endif // UNO_DOGLEGMETHOD_H