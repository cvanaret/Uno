// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNOPY_H
#define UNO_UNOPY_H

#include <functional>
#include <optional>
#include <vector>
#include "../UserModel.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/PointerWrapper.hpp"

namespace uno {
   //using Objective = std::function<double(const Vector<double>&)>;
   using Objective = double (*)(PointerWrapper<const Vector<double>>);

   using Constraints = std::function<void(const Vector<double>&, std::vector<double>&)>;

   //using ObjectiveGradient = std::function<void(const Vector<double>&, Vector<double>&)>;
   typedef void (*ObjectiveGradient)(PointerWrapper<const Vector<double>>, PointerWrapper<Vector<double>>);

   using Jacobian = std::function<void(const Vector<double>& /*x*/, double* /*jacobian_values*/)>;

   using Hessian = std::function<void(const Vector<double>& /*x*/, double objective_multiplier,
      const Vector<double>& /*multipliers*/, double* /*hessian_values*/)>;

   // TODO
   using JacobianOperator = void*;
   using JacobianTransposedOperator = void*;
   using HessianOperator = void*;

   using PythonUserModel = UserModel<Objective*, ObjectiveGradient, Constraints, Jacobian, JacobianOperator,
      JacobianTransposedOperator, Hessian, HessianOperator, std::optional<std::vector<double>>>;
} // namespace

#endif // UNO_UNOPY_H