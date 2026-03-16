// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LSR1HESSIAN_H
#define UNO_LSR1HESSIAN_H

#include "QuasiNewtonHessian.hpp"

namespace uno {
   // express the Hessian approximation at iteration k by a low-rank update:
   // Bk = B0 + U Uᵀ
   class LSR1Hessian: public QuasiNewtonHessian {
   public:
      LSR1Hessian(const Model& model, double objective_multiplier, const Options& options);
      ~LSR1Hessian() override = default;

      [[nodiscard]] bool is_positive_definite() const override;

      void initialize_statistics(Statistics& statistics) const override;
      void notify_accepted_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) override;
      void compute_hessian_vector_product(const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* result) override;

   protected:
      DenseMatrix<double> L; // lower triangular
      DenseMatrix<double> N;
      // Hessian representation: Bk = B0 + U P⁻¹ Uᵀ where B0 = delta I
      DenseMatrix<double> U;
      double delta{1.};
      bool hessian_recomputation_required{false};

      void recompute_hessian_representation();
   };
} // namespace

#endif // UNO_LSR1HESSIAN_H