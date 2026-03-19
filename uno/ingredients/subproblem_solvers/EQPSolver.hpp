// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EQPSOLVER_H
#define UNO_EQPSOLVER_H

#include <memory>
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "SubproblemSolver.hpp"

namespace uno {
   // forward declaration
   class Options;

   class EQPSolver: public SubproblemSolver {
   public:
      explicit EQPSolver(const Options& options);
      ~EQPSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius, const Vector<double>& initial_point,
         Direction& direction, Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   protected:
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
   };
} // namespace

#endif // UNO_EQPSOLVER_H