// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ZEROHESSIAN_H
#define UNO_ZEROHESSIAN_H

#include "HessianModel.hpp"

namespace uno {
   class ZeroHessian : public HessianModel {
   public:
      ZeroHessian() = default;

      void evaluate(const Model& model, const Vector<double>& primal_variables, double objective_multiplier,
         const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
   };
} // namespace

#endif // UNO_ZEROHESSIAN_H