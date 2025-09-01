// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PYTHONTYPES_H
#define UNO_PYTHONTYPES_H

#include <functional>
#include <vector>
#include "linear_algebra/Vector.hpp"
#include "tools/PointerWrapper.hpp"

namespace uno {
   // types of the model's functions (objective, constraints, gradients, Jacobian, Hessian)
   using objective_function_type = std::function<double(PointerWrapper<Vector<double>>)>;
   using constraint_functions_type = std::function<void(PointerWrapper<Vector<double>>, PointerWrapper<std::vector<double>>)>;
   using objective_gradient_type = std::function<void(PointerWrapper<Vector<double>>, PointerWrapper<Vector<double>>)>;
   using jacobian_type = std::function<void(PointerWrapper<Vector<double>> /*x*/,
      PointerWrapper<double> /*jacobian_values*/)>;
   using lagrangian_hessian_type = std::function<void(PointerWrapper<Vector<double>> /*x*/, double objective_multiplier,
      PointerWrapper<Vector<double>> /*y*/, PointerWrapper<double> /*hessian_values*/)>;
} // namespace

#endif // UNO_PYTHONTYPES_H