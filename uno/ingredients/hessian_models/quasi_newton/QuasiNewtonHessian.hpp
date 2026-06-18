// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QUASINEWTONHESSIAN_H
#define UNO_QUASINEWTONHESSIAN_H

#include "../HessianModel.hpp"
#include "linear_algebra/DenseMatrix.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declarations
   class Model;
   class Options;

   // express the Hessian approximation at iteration k by a low-rank update
   class QuasiNewtonHessian: public HessianModel {
   public:
      QuasiNewtonHessian(std::string_view name, const Model& model, double objective_multiplier, const Options& options);
      ~QuasiNewtonHessian() override = default;

      [[nodiscard]] bool has_hessian_operator() const override;
      [[nodiscard]] bool has_curvature() const override;

   protected:
      const Model& model;
      const double fixed_objective_multiplier;
      const size_t memory_size; // user defined
      // the (s, y) pairs are stored in CHRONOLOGICAL order in slots [0, number_entries_in_memory):
      // slot 0 holds the oldest pair, slot (number_entries_in_memory - 1) the newest. This ordering is mandatory:
      // the compact representation reproduces the BFGS matrix obtained by applying the updates in column order, and
      // since BFGS updates do not commute, a permuted (e.g. circular) slot order produces the wrong matrix.
      size_t number_entries_in_memory{0}; // 0 <= number_entries_in_memory <= memory_size
      // limited memory
      DenseMatrix<double> S;
      DenseMatrix<double> Y;
      Vector<double> current_lagrangian_gradient;
      Vector<double> trial_lagrangian_gradient;
      // most recent candidate pair (s, y), computed but not yet committed to the memory (see commit_memory_entry)
      Vector<double> latest_s;
      Vector<double> latest_y;
      double delta{1.};
      bool hessian_recomputation_required{false};
      const double delta_lower_bound;
      const double delta_upper_bound;

      // compute the candidate pair (s, y) into latest_s/latest_y, WITHOUT modifying the memory
      void compute_candidate_pair(const Iterate& current_iterate, const Iterate& trial_iterate, EvaluationCache& evaluation_cache);
      void compute_latest_s(const Iterate& current_iterate, const Iterate& trial_iterate);
      void compute_latest_y(const Iterate& current_iterate, const Iterate& trial_iterate, EvaluationCache& evaluation_cache);
      // append the candidate pair as the newest memory entry, physically dropping the oldest if the memory is full;
      // returns the slot of the new entry (always number_entries_in_memory - 1 afterwards)
      size_t commit_memory_entry();
      // physically drop the oldest entry (slot 0) and shift the remaining entries one slot toward the front; subclasses
      // extend this to shift their own per-entry cached data (e.g. the L matrix and the D diagonal)
      virtual void shift_memory_entries();
      // helpers for subclasses overriding shift_memory_entries: after the oldest entry is dropped, an entry indexed by
      // the pair (i, j) becomes (i-1, j-1), so a per-entry triangular cache is shifted in place by entry(i, j) =
      // entry(i+1, j+1) over the relevant triangle. include_diagonal selects whether the diagonal carries data.
      void shift_lower_triangle(DenseMatrix<double>& matrix, bool include_diagonal) const;
      void shift_upper_triangle(DenseMatrix<double>& matrix, bool include_diagonal) const;
      virtual void recompute_hessian_representation() = 0;
   };
} // namespace

#endif // UNO_QUASINEWTONHESSIAN_H
