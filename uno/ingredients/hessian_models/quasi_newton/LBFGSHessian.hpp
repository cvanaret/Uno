// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LBFGSHESSIAN_H
#define UNO_LBFGSHESSIAN_H

#include <vector>
#include "QuasiNewtonHessian.hpp"

namespace uno {
   // express the Hessian approximation at iteration k by a low-rank update:
   // Bk = B0 - U Uᵀ + V Vᵀ
   // where
   // B0 = δ I
   // V = Yk Dk^(-1/2)
   // M = δ Skᵀ Sk + Lk Dk⁻¹ Lkᵀ = J Jᵀ
   // U = (δ Sk + Yk Dk⁻¹ Lkᵀ) J⁻ᵀ
   class LBFGSHessian: public QuasiNewtonHessian {
   public:
      LBFGSHessian(const Model& model, double objective_multiplier, const Options& options);
      ~LBFGSHessian() override = default;

      [[nodiscard]] bool is_positive_definite() const override;

      void initialize_statistics(Statistics& statistics) const override;
      void notify_accepted_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) override;

      void compute_hessian_vector_product(const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* result) override;

      // functions that can be called by WoodburyEQPSolver
      [[nodiscard]] size_t get_correction_rank() const override;
      [[nodiscard]] VectorView<std::vector<double>> get_correction_column(size_t column_index) const override;
      [[nodiscard]] double get_correction_column_scaling(size_t column_index) const override;

   protected:
      DenseMatrix<double> L; // lower triangular
      std::vector<double> D; // diagonal
      std::vector<double> invsqrt_D; // diagonal
      DenseMatrix<double> L_invsqrt_D; // lower triangular
      DenseMatrix<double> M;
      // Hessian representation: Bk = B0 - U Uᵀ + V Vᵀ where B0 = δ I
      DenseMatrix<double> U;
      DenseMatrix<double> V;

      void update_D();
      void recompute_hessian_representation() override;
      [[nodiscard]] double compute_delta() const;
   };
} // namespace

#endif // UNO_LBFGSHESSIAN_H