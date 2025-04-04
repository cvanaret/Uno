// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LBFGSHESSIAN_H
#define UNO_LBFGSHESSIAN_H

#include "HessianModel.hpp"
#include "linear_algebra/DenseMatrix.hpp"

namespace uno {
   // forward declaration
   class Options;

   class LBFGSHessian : public HessianModel {
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
      void evaluate_hessian(Statistics& statistics, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* hessian_values) override;
      void compute_hessian_vector_product(const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* result) override;

   protected:
      const Model& model;
      const size_t memory_size;
      DenseMatrix<double> S_matrix{};
      DenseMatrix<double> Y_matrix{};
      DenseMatrix<double> M_matrix{};
   };
} // namespace

#endif // UNO_LBFGSHESSIAN_H