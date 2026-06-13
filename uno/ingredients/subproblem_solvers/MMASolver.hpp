// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MMASOLVER_H
#define UNO_MMASOLVER_H

#include <vector>
#include <memory>
#include "ingredients/subproblem_solvers/SubproblemSolver.hpp"
#include "ingredients/subproblem_solvers/SolverWorkspace.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class Options;

   class MMAWorkspace : public SolverWorkspace {
   public:
      MMAWorkspace() = default;
      ~MMAWorkspace() override = default;
      [[nodiscard]] double compute_hessian_quadratic_form(const Subproblem& /*subproblem*/, const Vector<double>& /*vector*/) const override {
         return 0.0; // MMA does not use quadratic Hessian models
      }
   };

   class MMASolver : public SubproblemSolver {
   public:
      MMASolver(size_t number_variables, size_t number_constraints, const Options& options);
      ~MMASolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& initial_point, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   private:
      size_t number_variables;
      size_t number_constraints;
      size_t iterations{0};

      // MMA history and parameters
      std::vector<double> x_old1;
      std::vector<double> x_old2;
      std::vector<double> lower_asymptotes;
      std::vector<double> upper_asymptotes;

      double asyinit;
      double asyincr;
      double asydecr;
      double external_move_limit;
      double internal_limit;
      size_t max_inner_iterations;

      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
      bool analysis_performed{false};

      void update_asymptotes(const Vector<double>& current_x, const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds);
   };
} // namespace

#endif // UNO_MMASOLVER_H
