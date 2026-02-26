// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LBFGSHESSIAN_H
#define UNO_LBFGSHESSIAN_H

#include <vector>
#include "../HessianModel.hpp"
#include "linear_algebra/DenseMatrix.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declarations
   class Model;
   class Options;

   // express the Hessian approximation at iteration k by a low-rank update:
   // Bk = B0 - U U^T + V V^T
   // where
   // B0 = delta_k I
   // V = Yk Dk^(-1/2)
   // U = (B0 Sk + Yk Dk^(-1) Lk^T) J^(-T)
   // J J^T = M = Sk^T B0 Sk + Lk Dk^(-1) Lk^T
   class LBFGSHessian: public HessianModel {
   public:
      LBFGSHessian(const Model& model, double objective_multiplier, const Options& options);
      ~LBFGSHessian() override = default;

      [[nodiscard]] bool has_hessian_operator() const override;
      [[nodiscard]] bool has_hessian_matrix() const override;
      [[nodiscard]] bool has_curvature() const override;
      [[nodiscard]] size_t number_nonzeros() const override;
      void compute_sparsity(int* row_indices, int* column_indices, int solver_indexing) const override;
      [[nodiscard]] bool is_positive_definite() const override;

      void initialize_statistics(Statistics& statistics) const override;
      void notify_accepted_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) override;
      void evaluate_hessian(Statistics& statistics, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* hessian_values) override;
      void compute_hessian_vector_product(const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* result) override;

   protected:
      const Model& model;
      const double fixed_objective_multiplier;
      const size_t memory_size; // user defined
      size_t number_entries_in_memory{0}; // 0 <= number_entries_in_memory <= memory_size
      size_t current_index{0}; // 0 <= current_index < memory_size
      // limited memory
      DenseMatrix<double> S;
      DenseMatrix<double> Y;
      DenseMatrix<double> L; // lower triangular
      std::vector<double> D; // diagonal
      std::vector<double> invsqrt_D; // diagonal
      DenseMatrix<double> L_invsqrt_D; // lower triangular
      DenseMatrix<double> M;
      // Hessian representation: Bk = B0 - U U^T + V V^T where B0 = delta I
      DenseMatrix<double> U;
      DenseMatrix<double> V;
      Vector<double> current_lagrangian_gradient;
      Vector<double> trial_lagrangian_gradient;
      double delta{1.};
      bool hessian_recomputation_required{false};

      void update_S(const Iterate& current_iterate, const Iterate& trial_iterate);
      void update_Y(const Iterate& current_iterate, const Iterate& trial_iterate, EvaluationCache& evaluation_cache);
      void update_D();
      void recompute_hessian_representation();
      [[nodiscard]] double compute_delta() const;
   };
} // namespace

#endif // UNO_LBFGSHESSIAN_H