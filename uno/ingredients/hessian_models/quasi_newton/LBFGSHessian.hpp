// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LBFGSHESSIAN_H
#define UNO_LBFGSHESSIAN_H

#include <vector>
#include "../HessianModel.hpp"
#include "linear_algebra/DenseMatrix.hpp"

namespace uno {
   class LBFGSHessian: public HessianModel {
   public:
      LBFGSHessian(size_t number_variables, size_t memory_size);
      ~LBFGSHessian() override = default;

      void initialize_statistics(Statistics& statistics, const Options& options) const override;
      void notify_accepted_iterate(const Model& model, const Iterate& current_iterate, const Iterate& trial_iterate) override;
      void evaluate_hessian(Statistics& statistics, const Model& model, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
      void compute_hessian_vector_product(const Model& model, const Vector<double>& vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, Vector<double>& result) override;

   protected:
      const size_t number_variables;
      const size_t memory_size; // user defined
      size_t used_memory_size{0}; // 0 <= used_memory_size <= memory_size
      size_t current_memory_slot; // 0 <= current_available_slot < memory_size
      // memory
      DenseMatrix<double> S_matrix;
      DenseMatrix<double> Y_matrix;
      // Hessian representation
      bool hessian_recomputation_required{false};
      DenseMatrix<double> L_matrix;
      std::vector<double> D_matrix; // diagonal
      DenseMatrix<double> M_matrix;

      void update_memory(const Model& model, const Iterate& current_iterate, const Iterate& trial_iterate);
      void recompute_hessian_representation();
      static void perform_high_rank_update(DenseMatrix<double>& matrix, size_t matrix_dimension, size_t matrix_leading_dimension,
         DenseMatrix<double>& high_rank_correction, size_t correction_rank, size_t correction_leading_dimension);
      static void perform_high_rank_update_transpose(DenseMatrix<double>& matrix, size_t matrix_dimension, size_t matrix_leading_dimension,
         DenseMatrix<double>& high_rank_correction, size_t correction_rank, size_t correction_leading_dimension);
      static void compute_cholesky_factors(DenseMatrix<double>& matrix, size_t dimension, size_t leading_dimension);
   };
} // namespace

#endif // UNO_LBFGSHESSIAN_H