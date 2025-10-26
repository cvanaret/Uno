// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LBFGSHESSIAN_H
#define UNO_LBFGSHESSIAN_H

#include <vector>
#include "../HessianModel.hpp"
#include "linear_algebra/DenseMatrix.hpp"

namespace uno {
   // forward declaration
   class Options;

   class LBFGSHessian: public HessianModel {
   public:
      LBFGSHessian(const Model& model, const Options& options);
      ~LBFGSHessian() override = default;

      [[nodiscard]] bool has_hessian_operator() const override;
      [[nodiscard]] bool has_hessian_matrix() const override;
      [[nodiscard]] bool has_curvature() const override;
      [[nodiscard]] size_t number_nonzeros() const override;
      void compute_sparsity(int* row_indices, int* column_indices, int solver_indexing) const override;
      [[nodiscard]] bool is_positive_definite() const override;

      void initialize(const Model& model) override;
      void initialize_statistics(Statistics& statistics, const Options& options) const override;
      void notify_accepted_iterate(Iterate& current_iterate, Iterate& trial_iterate) override;
      void evaluate_hessian(Statistics& statistics, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* hessian_values) override;
      void compute_hessian_vector_product(const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* result) override;

   protected:
      const Model& model;
      const size_t memory_size; // user defined
      size_t current_memory_size{0}; // 0 <= used_memory_size <= memory_size
      size_t current_available_slot{0}; // 0 <= current_available_slot < memory_size
      // memory
      DenseMatrix<double> S_matrix;
      DenseMatrix<double> Y_matrix;
      // Hessian representation
      bool hessian_recomputation_required{false};
      DenseMatrix<double> L_matrix;
      std::vector<double> D_matrix; // diagonal
      DenseMatrix<double> M_matrix;

      void update_memory(const Model& model, Iterate& current_iterate, Iterate& trial_iterate);
      void recompute_hessian_representation();
      static void perform_high_rank_update(DenseMatrix<double>& matrix, size_t matrix_dimension, size_t matrix_leading_dimension,
         DenseMatrix<double>& high_rank_correction, size_t correction_rank, size_t correction_leading_dimension);
      static void perform_high_rank_update_transpose(DenseMatrix<double>& matrix, size_t matrix_dimension, size_t matrix_leading_dimension,
         DenseMatrix<double>& high_rank_correction, size_t correction_rank, size_t correction_leading_dimension);
      static void compute_cholesky_factors(DenseMatrix<double>& matrix, size_t dimension, size_t leading_dimension);
   };
} // namespace

#endif // UNO_LBFGSHESSIAN_H