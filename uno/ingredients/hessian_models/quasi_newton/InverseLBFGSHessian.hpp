// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INVERSELBFGSHESSIAN_H
#define UNO_INVERSELBFGSHESSIAN_H

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
   class InverseLBFGSHessian: public QuasiNewtonHessian {
   public:
      InverseLBFGSHessian(const Model& model, double objective_multiplier, const Options& options);
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
      DenseMatrix<double> L; // lower triangular
      std::vector<double> D; // diagonal
      std::vector<double> invsqrt_D; // diagonal
      DenseMatrix<double> L_invsqrt_D; // lower triangular
      DenseMatrix<double> M;
      // Hessian representation: Bk = B0 - U Uᵀ + V Vᵀ where B0 = δ I
      DenseMatrix<double> U;
      DenseMatrix<double> V;
      const double delta_upper_bound;

      void update_D();
      void recompute_hessian_representation() override;
      [[nodiscard]] double compute_delta() const;
   };
} // namespace

#endif // UNO_INVERSELBFGSHESSIAN_H