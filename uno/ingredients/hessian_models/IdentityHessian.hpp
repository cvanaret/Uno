// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_IDENTITYHESSIAN_H
#define UNO_IDENTITYHESSIAN_H

#include "HessianModel.hpp"

namespace uno {
   class IdentityHessian : public HessianModel {
   public:
      IdentityHessian() = default;

      void initialize(const Model& model) override;
      [[nodiscard]] size_t number_nonzeros(const Model& model) const override;
      [[nodiscard]] bool is_positive_definite() const override;
      void evaluate_hessian(Statistics& statistics, const Model& model, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
      void compute_hessian_vector_product(const Model& model, const Vector<double>& vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, Vector<double>& result) override;
      [[nodiscard]] std::string get_name() const override;
   };
} // namespace

#endif // UNO_IDENTITYHESSIAN_H