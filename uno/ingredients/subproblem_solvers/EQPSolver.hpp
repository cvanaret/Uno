// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EQPSOLVER_H
#define UNO_EQPSOLVER_H

#include <memory>
#include <string>
#include "LowRankCorrection.hpp"
#include "SubproblemSolver.hpp"
#include "optimization/Direction.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class NewtonSystem;
   class Options;
   class Statistics;

   // iterative refinement on the uncondensed Newton system (PDFullSpaceSolver::Solve in IPOPT)
   struct IterativeRefinementParameters {
      size_t min_steps; // minimum number of refinement steps
      size_t max_steps; // maximum number of refinement steps (0: iterative refinement disabled)
      double residual_ratio_max; // refinement stops below this residual ratio
      double residual_ratio_singular; // above this residual ratio after failed refinement, the matrix is deemed singular
      double improvement_factor; // minimal required reduction of the residual ratio per refinement step
   };

   // LowRankCorrection = NoLowRankCorrection: explicit Hessian in the augmented matrix
   // LowRankCorrection = WoodburyCorrection: diagonal + low-rank quasi-Newton Hessian (L-BFGS, L-SR1)
   template <typename LowRankCorrection>
   class EQPSolver: public SubproblemSolver {
   public:
      explicit EQPSolver(const Options& options, LowRankCorrection correction = {});
      ~EQPSolver() override; // defined in EQPSolver.cpp (incomplete linear solver type)

      void initialize_memory(const Subproblem& subproblem) override;
      void initialize_statistics(Statistics& statistics) override;
      void compute_least_squares_multipliers(const Subproblem& subproblem, Iterate& iterate, Evaluations& evaluations,
         double multipliers_threshold) override;

      [[nodiscard]] const Direction& solve(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& initial_point, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] bool has_second_order_corrections() const override;
      void initialize_second_order_corrections(const Subproblem& subproblem, const Iterate& current_iterate,
         const Iterate& trial_iterate, Evaluations& current_evaluations, Evaluations& trial_evaluations) override;
      const Direction& compute_second_order_correction(Statistics& statistics, const Subproblem& subproblem,
         const Iterate& current_iterate) override;
      void update_second_order_corrections(const Subproblem& subproblem, const Iterate& trial_iterate,
         Evaluations& trial_evaluations) override;

      [[nodiscard]] const SolverWorkspace& get_workspace() const override;

      [[nodiscard]] const NewtonSystem& get_newton_system() const { return *this->newton_system; }

      // diagnostics of the last solve (exposed for testing). The residual ratio is NaN when refinement is disabled
      [[nodiscard]] double get_residual_ratio() const { return this->residual_ratio; }
      [[nodiscard]] size_t get_number_refinement_steps() const { return this->number_refinement_steps; }
      [[nodiscard]] bool is_matrix_deemed_singular() const { return this->matrix_deemed_singular; }

   protected:
      Direction direction;
      std::unique_ptr<NewtonSystem> newton_system;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
      LowRankCorrection low_rank_correction;
      bool analysis_performed{false};
      Vector<double> zero_multipliers{};
      // the least-squares system overwrote the augmented matrix and its factorization: the next direction computation
      // must assemble and factorize the augmented matrix again
      bool factorization_holds_least_squares_system{false};
      // regularization of the (2, 2) block of the least-squares system when the Jacobian is rank deficient
      static constexpr double least_squares_dual_regularization{1e-8};

      // iterative refinement
      IterativeRefinementParameters iterative_refinement; // not const: adjustable by derived classes (tests)
      Vector<double> uncondensed_rhs{};
      Vector<double> uncondensed_solution{};
      Vector<double> uncondensed_residual{};
      Vector<double> uncondensed_correction{};
      Vector<double> condensed_rhs_buffer{};
      Vector<double> condensed_solution_buffer{};
      // whether the pivoting tolerance was increased for the current matrix (IPOPT's augsys_improved_). The tolerance
      // itself is never decreased
      bool pivoting_tolerance_increased_for_current_matrix{false};
      double residual_ratio{0.};
      size_t number_refinement_steps{0};
      bool matrix_deemed_singular{false}; // in the last solve, the matrix was deemed singular (IPOPT's "pretend singular")
      // per-iteration summary (IPOPT-style): refinement steps and flags (scal: scaling switched on, piv: pivoting tolerance
      // increased, sing: matrix deemed singular, inac: inaccurate solution accepted). Shown in the extended "IR" column
      size_t iteration_refinement_steps{0};
      std::vector<std::string_view> iteration_refinement_flags{};

      void solve_with_iterative_refinement(Statistics& statistics, const Subproblem& subproblem,
         bool allow_refactorization);
      [[nodiscard]] bool refine_uncondensed_solution();
      [[nodiscard]] double compute_uncondensed_residual_ratio();
      void compute_iterative_refinement_correction();
      void report_iterative_refinement(Statistics& statistics);

      bool SOC_initialized{false};
      Vector<double> constraints_SOC;
      Vector<double> constraints_buffer_SOC;
      Direction direction_SOC;
   };

   // explicitly instantiated in EQPSolver.cpp
   extern template class EQPSolver<NoLowRankCorrection>;
   extern template class EQPSolver<WoodburyCorrection>;
} // namespace

#endif // UNO_EQPSOLVER_H