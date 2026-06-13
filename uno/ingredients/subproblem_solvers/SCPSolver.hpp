// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCPSOLVER_H
#define UNO_SCPSOLVER_H

#include <vector>
#include <memory>
#include "ingredients/subproblem_solvers/SubproblemSolver.hpp"
#include "ingredients/subproblem_solvers/SolverWorkspace.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"

namespace uno {
   class Options;

   class SCPWorkspace : public SolverWorkspace {
   public:
      SCPWorkspace() = default;
      ~SCPWorkspace() override = default;
      [[nodiscard]] double compute_hessian_quadratic_form(const Subproblem& /*subproblem*/, const Vector<double>& /*vector*/) const override {
         return 0.0; // SCP approximations have strictly diagonal/separable Hessians, mapped into the linear system natively
      }
   };

   class SCPSolver : public SubproblemSolver {
   public:
      SCPSolver(size_t number_variables, size_t number_constraints, const Options& options);
      ~SCPSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& initial_point, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   protected:
      size_t number_variables;
      size_t number_constraints;
      size_t iterations{0};

      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
      bool analysis_performed{false};

      // Pure virtual method for derived classes to inject their specific convex approximation diagonals
      virtual void compute_diagonal_hessian(const Vector<double>& initial_point, Evaluations& current_evaluations, std::vector<double>& Dx) = 0;
   };
} // namespace

#endif // UNO_SCPSOLVER_H
