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
   class Options;
   class Statistics;

   // iterative refinement on the uncondensed Newton system (PDFullSpaceSolver::Solve in IPOPT)
   struct IterativeRefinementParameters {
      size_t min_steps; // minimum number of refinement steps
      size_t max_steps; // maximum number of refinement steps (0: iterative refinement disabled)
      double residual_ratio_max; // refinement stops below this residual ratio
      double residual_ratio_singular; // above this residual ratio after failed refinement, the matrix is deemed singular
      double improvement_factor; // minimal required reduction of the residual ratio per refinement step
      // undo the last pivoting tolerance increase after this many consecutive iterations whose (unrefined) solution
      // was already accurate (0: never, as in IPOPT). Doubled each time the tolerance must be increased again soon
      // after a decrease (oscillation backoff)
      size_t pivoting_tolerance_decrease_after;
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
      const Direction& compute_second_order_correction(const Subproblem& subproblem, const Iterate& current_iterate) override;
      void update_second_order_corrections(const Subproblem& subproblem, const Iterate& trial_iterate,
         Evaluations& trial_evaluations) override;

      [[nodiscard]] const SolverWorkspace& get_workspace() const override;

      // diagnostics of the last solve (exposed for testing). The residual ratio is NaN when refinement is disabled
      [[nodiscard]] double get_residual_ratio() const { return this->residual_ratio; }
      [[nodiscard]] size_t get_number_refinement_steps() const { return this->number_refinement_steps; }
      [[nodiscard]] bool is_matrix_deemed_singular() const { return this->matrix_deemed_singular; }
      // number of pivoting tolerance increases since the creation of the solver
      [[nodiscard]] size_t get_number_pivoting_tolerance_increases() const { return this->number_pivoting_tolerance_increases; }
      [[nodiscard]] size_t get_number_pivoting_tolerance_decreases() const { return this->number_pivoting_tolerance_decreases; }

   protected:
      Direction direction;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
      LowRankCorrection correction;
      bool analysis_performed{false};

      // iterative refinement
      IterativeRefinementParameters iterative_refinement; // not const: adjustable by derived classes (tests)
      Vector<double> uncondensed_rhs{};
      Vector<double> uncondensed_solution{};
      Vector<double> uncondensed_residual{};
      Vector<double> uncondensed_correction{};
      Vector<double> augmented_rhs_buffer{};
      Vector<double> augmented_solution_buffer{};
      // whether the pivoting tolerance was increased for the current matrix (IPOPT's augsys_improved_)
      bool pivoting_tolerance_increased_for_current_matrix{false};
      size_t number_consecutive_accurate_iterations{0};
      size_t pivoting_tolerance_decrease_after;
      size_t number_iterations_since_decrease{0};
      bool pivoting_tolerance_decreased{false};
      bool current_iteration_accurate{false}; // no previous iteration
      double residual_ratio{0.};
      size_t number_refinement_steps{0};
      bool matrix_deemed_singular{false};
      size_t number_pivoting_tolerance_increases{0};
      size_t number_pivoting_tolerance_decreases{0};
      // per-iteration summary (IPOPT-style): refinement steps, q: pivoting tolerance increased, s: matrix deemed singular,
      // d: pivoting tolerance decreased, S: inaccurate solution accepted. Shown in the extended "IR" statistics column
      Statistics* statistics{nullptr};
      size_t iteration_refinement_steps{0};
      std::string iteration_refinement_flags{};
      void report_iterative_refinement();

      void solve_with_iterative_refinement(Statistics* statistics, const Subproblem& subproblem, const Iterate& current_iterate);
      [[nodiscard]] double compute_uncondensed_residual_ratio(const Subproblem& subproblem, const Iterate& current_iterate);
      void compute_iterative_refinement_correction(const Subproblem& subproblem, const Iterate& current_iterate);

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