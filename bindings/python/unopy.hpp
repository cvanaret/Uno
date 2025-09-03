// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNOPY_H
#define UNO_UNOPY_H

#include <functional>
#include "../UserModel.hpp"
#include "tools/PointerWrapper.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   using Objective = std::function<double(PointerWrapper<Vector<double>>)>;

   using Constraints = std::function<void(PointerWrapper<Vector<double>>, PointerWrapper<std::vector<double>>)>;

   using ObjectiveGradient = std::function<void(PointerWrapper<Vector<double>>, PointerWrapper<Vector<double>>)>;

   using Jacobian = std::function<void(PointerWrapper<Vector<double>> /*x*/,
      PointerWrapper<double> /*jacobian_values*/)>;

   using Hessian = std::function<void(PointerWrapper<Vector<double>> /*x*/, double objective_multiplier,
      PointerWrapper<Vector<double>> /*y*/, PointerWrapper<double> /*hessian_values*/)>;

   // TODO
   using JacobianOperator = void*;
   using JacobianTransposedOperator = void*;
   using HessianOperator = void*;

   using PythonUserModel = UserModel<Objective, ObjectiveGradient, Constraints, Jacobian, JacobianOperator,
      JacobianTransposedOperator, Hessian, HessianOperator>;
} // namespace

#endif // UNO_UNOPY_H