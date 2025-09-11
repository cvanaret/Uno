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
   using Objective = std::function<int(const Vector<double>&, PointerWrapper<double>)>;

   using Constraints = std::function<int(const Vector<double>&, PointerWrapper<double>)>;

   using ObjectiveGradient = std::function<int(const Vector<double>&, PointerWrapper<double>)>;

   using Jacobian = std::function<int(const Vector<double>&, PointerWrapper<double>)>;

   using Hessian = std::function<int(const Vector<double>&, double objective_multiplier, const Vector<double>&,
      PointerWrapper<double>)>;

   using JacobianOperator = void*;
   using JacobianTransposedOperator = void*;

   using HessianOperator = std::function<int(PointerWrapper<const double>, bool, double objective_multiplier,
      const Vector<double>&, PointerWrapper<const double>, PointerWrapper<double>)>;

   using PythonUserModel = UserModel<std::optional<Objective>, std::optional<ObjectiveGradient>, std::optional<Constraints>,
      std::optional<Jacobian>, JacobianOperator, JacobianTransposedOperator, std::optional<Hessian>, std::optional<HessianOperator>,
      std::optional<std::vector<double>>>;
} // namespace

#endif // UNO_UNOPY_H