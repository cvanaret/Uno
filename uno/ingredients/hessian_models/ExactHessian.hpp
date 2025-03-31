// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EXACTHESSIAN_H
#define UNO_EXACTHESSIAN_H

#include "HessianModel.hpp"

namespace uno {
   class ExactHessian : public HessianModel {
   public:
      ExactHessian() = default;

      [[nodiscard]] size_t compute_number_hessian_nonzeros(const Model& model) const override;
      void evaluate(const Model& model, const Vector<double>& primal_variables, double objective_multiplier,
         const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
   };
} // namespace

#endif // UNO_EXACTHESSIAN_H