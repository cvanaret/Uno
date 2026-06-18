// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INVERSELBFGSHESSIAN_H
#define UNO_INVERSELBFGSHESSIAN_H

#include "../QuasiNewtonHessian.hpp"

namespace uno {
   // express the inverse Hessian approximation Hk = Bk⁻¹ at iteration k by a low-rank update (compact inverse BFGS
   // representation), using the upper triangular matrix R with R(i, j) = sᵢᵀ yⱼ for i ≤ j and H0 = δ⁻¹ I.
   // Sk, Yk and R are expressed in chronological order (slot 0 oldest, last slot newest); when the memory is full, the
   // oldest entry is physically shifted out so that the slot order always matches the chronological order. This ordering
   // is required: R is the upper triangular part of Skᵀ Yk in chronological order, and the BFGS updates do not commute.
   class InverseLBFGSHessian: public QuasiNewtonHessian {
   public:
      InverseLBFGSHessian(const Model& model, const Options& options);
      ~InverseLBFGSHessian() override = default;

      [[nodiscard]] bool has_hessian_matrix() const override;
      [[nodiscard]] size_t number_nonzeros() const override;
      void compute_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing) const override;
      [[nodiscard]] bool is_positive_definite() const override;

      void initialize_statistics(Statistics& statistics) const override;
      void notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) override;

      void evaluate_hessian(Statistics& statistics, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* hessian_values) override;
      void compute_hessian_vector_product(const double* x, const double* vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, double* result) override;

      // function that can be called by NewtonSolver
      void compute_inverse_hessian_vector_product(const double* x, const double* vector, double* result);

   protected:
      DenseMatrix<double> R; // upper triangular, R(i, j) = sᵢᵀ yⱼ for i ≤ j (diagonal included)
      // temporaries
      Vector<double> STv; // Sᵀv
      Vector<double> YTv; // Yᵀv

      void shift_memory_entries() override;
      void recompute_hessian_representation() override;
      [[nodiscard]] double compute_delta() const;
   };
} // namespace

#endif // UNO_INVERSELBFGSHESSIAN_H
