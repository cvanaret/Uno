// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EXACTHESSIAN_H
#define UNO_EXACTHESSIAN_H

#include "HessianModel.hpp"

namespace uno {
   class ExactHessian : public HessianModel {
   public:
      ExactHessian() = default;
      ~ExactHessian() override = default;

      [[nodiscard]] bool has_hessian_operator(const Model& model) const override;
      [[nodiscard]] bool has_hessian_matrix(const Model& model) const override;
      [[nodiscard]] bool has_curvature(const Model& model) const override;
      [[nodiscard]] size_t number_nonzeros(const Model& model) const override;
      void compute_sparsity(const Model& model, int* row_indices, int* column_indices, int solver_indexing) const override;
      [[nodiscard]] bool is_positive_definite() const override;

      void initialize(const Model& model) override;
      void initialize_statistics(Statistics& statistics, const Options& options) const override;
      void evaluate_hessian(Statistics& statistics, const Model& model, const Vector<double>& primal_variables, double objective_multiplier,
         const Vector<double>& constraint_multipliers, double* hessian_values) override;
      void compute_hessian_vector_product(const Model& model, const double* x, const double* vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, double* result) override;
      [[nodiscard]] std::string get_name() const override;
   };
} // namespace

#endif // UNO_EXACTHESSIAN_H