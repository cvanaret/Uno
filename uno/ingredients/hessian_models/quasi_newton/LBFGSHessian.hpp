// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LBFGSHESSIAN_H
#define UNO_LBFGSHESSIAN_H

#include <vector>
#include "../HessianModel.hpp"
#include "linear_algebra/DenseMatrix.hpp"

namespace uno {
   // express the Hessian approximation at iteration k by a low-rank update:
   // Bk = B0 - U U^T + V V^T
   // where
   // B0 = delta_k I
   // V = Yk Dk^(-1/2)
   // U = (B0 Sk + Yk Dk^(-1) Lk^T) J^(-T)
   // J J^T = M = Sk^T B0 Sk + Lk Dk^(-1) Lk^T
   class LBFGSHessian: public HessianModel {
   public:
      explicit LBFGSHessian(const Options& options);
      ~LBFGSHessian() override = default;

      [[nodiscard]] bool has_implicit_representation() const override;
      [[nodiscard]] bool has_explicit_representation() const override;
      [[nodiscard]] bool has_curvature(const Model& model) const override;
      [[nodiscard]] size_t number_nonzeros(const Model& model) const override;
      [[nodiscard]] bool is_positive_definite() const override;

      void initialize(const Model& model) override;
      void initialize_statistics(Statistics& statistics, const Options& options) const override;
      void notify_accepted_iterate(const Model& model, Iterate& current_iterate, Iterate& trial_iterate) override;
      void evaluate_hessian(Statistics& statistics, const Model& model, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
      void compute_hessian_vector_product(const Model& model, const double* vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, double* result) override;
      [[nodiscard]] std::string get_name() const override;

   protected:
      size_t dimension{};
      const size_t memory_size; // user defined
      size_t number_entries_in_memory{0}; // 0 <= used_memory_size <= memory_size
      size_t current_memory_slot{0}; // 0 <= current_available_slot < memory_size
      // memory
      DenseMatrix<double> S_matrix;
      DenseMatrix<double> Y_matrix;
      // Hessian representation: Bk = B0 - U U^T + V V^T
      bool hessian_recomputation_required{false};
      DenseMatrix<double, MatrixShape::LOWER_TRIANGULAR> L_matrix;
      std::vector<double> D_matrix; // diagonal
      DenseMatrix<double> M_matrix;
      DenseMatrix<double> U_matrix;
      DenseMatrix<double> V_matrix;
      double initial_identity_multiple{1.}; // referred to as delta in Numerical optimization

      void update_memory(const Model& model, Iterate& current_iterate, Iterate& trial_iterate);
      void update_Y_matrix(const Model& model, Iterate& current_iterate, Iterate& trial_iterate);
      void update_S_matrix(const Iterate& current_iterate, const Iterate& trial_iterate);
      void update_D_matrix();
      void recompute_hessian_representation();
      [[nodiscard]] double compute_initial_identity_factor() const;

      static void perform_high_rank_update(DenseMatrix<double>& matrix, size_t matrix_dimension, size_t matrix_leading_dimension,
         DenseMatrix<double, MatrixShape::LOWER_TRIANGULAR>& high_rank_correction, size_t correction_rank, size_t correction_leading_dimension, double alpha, double beta);
      static void perform_high_rank_update_transpose(DenseMatrix<double>& matrix, size_t matrix_dimension, size_t matrix_leading_dimension,
         DenseMatrix<double>& high_rank_correction, size_t correction_rank, size_t correction_leading_dimension, double alpha, double beta);
      static void compute_cholesky_factors(DenseMatrix<double>& matrix, size_t dimension, size_t leading_dimension);
      static void lower_triangular_back_solve(DenseMatrix<double>& L_matrix, size_t matrix_leading_dimension, DenseMatrix<double>& rhs,
         size_t rhs_dimension, size_t leading_dimension_rhs, bool transpose);
   };
} // namespace

#endif // UNO_LBFGSHESSIAN_H