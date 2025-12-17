// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_IDENTITYHESSIAN_H
#define UNO_IDENTITYHESSIAN_H

#include "HessianModel.hpp"

namespace uno {
   class IdentityHessian : public HessianModel {
   public:
      explicit IdentityHessian(size_t number_variables);

      [[nodiscard]] bool has_hessian_operator() const override;
      [[nodiscard]] bool has_hessian_matrix() const override;
      [[nodiscard]] bool has_curvature() const override;
      [[nodiscard]] size_t number_nonzeros() const override;
      void compute_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing) const override;
      [[nodiscard]] bool is_positive_definite() const override;

      void evaluate_hessian(Statistics& statistics, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* hessian_values) override;
      void compute_hessian_vector_product(const double* x, const double* vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, double* result) override;

   protected:
      const size_t number_variables;
   };
} // namespace

#endif // UNO_IDENTITYHESSIAN_H