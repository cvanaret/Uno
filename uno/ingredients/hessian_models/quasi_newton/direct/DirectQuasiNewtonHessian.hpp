// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DIRECTQUASINEWTONHESSIAN_H
#define UNO_DIRECTQUASINEWTONHESSIAN_H

#include "../QuasiNewtonHessian.hpp"
#include "linear_algebra/DenseMatrix.hpp"

namespace uno {
   class DirectQuasiNewtonHessian: public QuasiNewtonHessian {
   public:
      DirectQuasiNewtonHessian(const std::string_view name, const Model& model, double objective_multiplier, const Options& options);
      ~DirectQuasiNewtonHessian() override = default;

      [[nodiscard]] bool has_hessian_operator() const override;
      [[nodiscard]] bool has_hessian_matrix() const override;
      [[nodiscard]] bool has_curvature() const override;
      [[nodiscard]] size_t number_nonzeros() const override;
      void compute_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing) const override;

      void evaluate_hessian(Statistics& statistics, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* hessian_values) override;

      // functions that can be called by WoodburyEQPSolver
      [[nodiscard]] virtual size_t get_correction_rank() const = 0;
      [[nodiscard]] virtual VectorView<std::vector<double>> get_correction_column(size_t column_index) const = 0;
      [[nodiscard]] virtual double get_correction_column_scaling(size_t column_index) const = 0;
   };
} // namespace

#endif // UNO_DIRECTQUASINEWTONHESSIAN_H