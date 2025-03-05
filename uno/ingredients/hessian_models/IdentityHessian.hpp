// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_IDENTITYHESSIAN_H
#define UNO_IDENTITYHESSIAN_H

#include "HessianModel.hpp"

namespace uno {
   class IdentityHessian : public HessianModel {
   public:
      IdentityHessian() = default;

      void evaluate(const Model& model, const Vector<double>& primal_variables, double objective_multiplier,
         const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
   };
} // namespace

#endif // UNO_IDENTITYHESSIAN_H