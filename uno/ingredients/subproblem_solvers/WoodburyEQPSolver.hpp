// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_WOODBURYEQPSOLVER_H
#define UNO_WOODBURYEQPSOLVER_H

#include <memory>
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "SubproblemSolver.hpp"

namespace uno {
   // forward declarations
   class LBFGSHessian;
   class Options;

   class WoodburyEQPSolver: public SubproblemSolver {
   public:
      WoodburyEQPSolver(LBFGSHessian& hessian_model, const Options& options);
      ~WoodburyEQPSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius, const Vector<double>& initial_point,
         Direction& direction, Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   protected:
      LBFGSHessian& hessian_model;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
   };
} // namespace

#endif // UNO_WOODBURYEQPSOLVER_H