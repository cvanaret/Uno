// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCALINGSTRATEGY_H
#define UNO_SCALINGSTRATEGY_H

#include <cstddef>
#include <optional>
#include <vector>
#include "SymmetricMatrixScaling.hpp"

namespace uno {
   // forward declaration
   class COOLinearSystem;

   struct LinearSystemScalingParameters {
      LinearSystemScalingMethod method{LinearSystemScalingMethod::RUIZ};
      // true: the scaling is switched on only when the factorizations seem inaccurate (IPOPT's linear_scaling_on_demand),
      // then kept until the end. false: all the matrices are scaled
      bool on_demand{true};
      // the scaling is switched on (from the next matrix on) when the first factorization of a matrix delays more than
      // this fraction of the pivots. >= 1: disabled (default: the delayed pivots of KKT matrices are mostly structural)
      double delayed_pivots_threshold{1.};
   };

   // symmetric diagonal scaling S K S of the matrix K of a direct linear solver (DirectSymmetricIndefiniteLinearSolver
   // applies it around the factorization and the solves: S K S y = S b, x = S y). The callers keep reading and writing the
   // unscaled K (e.g. for the residuals of iterative refinement): the scaled values only exist during the factorization.
   // Cost control:
   // - on demand, nothing is done (zero overhead) until the scaling is switched on, either by the caller (activate, e.g.
   //   iterative refinement failed) or because the factorization delayed too many pivots;
   // - the scaling factors are computed at most once per matrix: at the first factorization after begin_new_matrix(). The
   //   refactorizations of the same matrix (inertia correction retries, pivoting tolerance increase) reuse them, and so do
   //   the solves with an existing factorization (second-order corrections, iterative refinement)
   class ScalingStrategy {
   public:
      explicit ScalingStrategy(const LinearSystemScalingParameters& parameters);

      // the next factorization is that of a new matrix (as opposed to a refactorization of the current matrix with a
      // different regularization or pivoting tolerance)
      void begin_new_matrix();
      // switch on the scaling from the next factorization on. Returns false if already on
      [[nodiscard]] bool activate();
      [[nodiscard]] bool is_active() const { return this->active; }

      // factorization: scale the matrix values in place if the scaling is on (returns whether they were scaled), then
      // restore them after the factorization (also when it throws)
      [[nodiscard]] bool scale_matrix(COOLinearSystem& linear_system);
      void unscale_matrix(COOLinearSystem& linear_system) const;
      // after the factorization: monitor the delayed pivots of the first unscaled factorization of a matrix
      void notify_factorization(std::optional<size_t> number_delayed_pivots, size_t dimension);

      // solves with a scaled factorization: S b (into an internal buffer) and x = S y (in place). Column-major blocks
      [[nodiscard]] bool is_last_factorization_scaled() const { return this->last_factorization_scaled; }
      [[nodiscard]] const double* scale_rhs(const double* rhs, size_t dimension, size_t number_of_rhs);
      void scale_solution(double* solution, size_t dimension, size_t number_of_rhs) const;

      // diagnostics (exposed for testing)
      [[nodiscard]] size_t get_number_scaling_computations() const { return this->number_scaling_computations; }
      [[nodiscard]] const std::vector<double>& get_scaling() const { return this->scaling; }

   protected:
      const LinearSystemScalingParameters parameters;
      SymmetricMatrixScaling scaling_method{};
      bool active;
      bool pending_activation{false}; // the scaling will be switched on at the next matrix
      bool new_matrix{true}; // the next factorization is the first factorization of a new matrix
      bool first_factorization_of_matrix{false}; // the current factorization is the first of its matrix
      bool scaling_up_to_date{false}; // the scaling was computed for the current matrix
      bool last_factorization_scaled{false};
      size_t number_scaling_computations{0};
      std::vector<double> scaling{};
      std::vector<double> entry_factors{};
      std::vector<double> unscaled_values{};
      std::vector<double> scaled_rhs{};

      void compute_scaling(const COOLinearSystem& linear_system);
   };
} // namespace

#endif // UNO_SCALINGSTRATEGY_H
