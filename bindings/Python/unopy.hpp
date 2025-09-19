// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNOPY_H
#define UNO_UNOPY_H

#include <functional>
#include <optional>
#include <vector>
#include <pybind11/pybind11.h>
#include "../UserModel.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/PointerWrapper.hpp"

namespace py = pybind11;

namespace uno {
   using Objective = std::function<int(int32_t, const Vector<double>&, PointerWrapper<double>, const py::object&)>;

   using Constraints = std::function<int(int32_t, int32_t, const Vector<double>&, PointerWrapper<double>, const py::object&)>;

   using ObjectiveGradient = std::function<int(int32_t, const Vector<double>&, PointerWrapper<double>, const py::object&)>;

   using Jacobian = std::function<int(int32_t, int32_t, const Vector<double>&, PointerWrapper<double>, const py::object&)>;

   using Hessian = std::function<int(int32_t, int32_t, int32_t, const Vector<double>&, double, const Vector<double>&,
      PointerWrapper<double>, const py::object&)>;

   using JacobianOperator = std::function<int(int32_t, int32_t, PointerWrapper<const double>, bool, PointerWrapper<const double>,
      PointerWrapper<double>, const py::object&)>;

   using JacobianTransposedOperator = std::function<int(int32_t, int32_t, PointerWrapper<const double>, bool,
      PointerWrapper<const double>, PointerWrapper<double>, const py::object&)>;

   using HessianOperator = std::function<int(int32_t, int32_t, PointerWrapper<const double>, bool, double,
      const Vector<double>&, PointerWrapper<const double>, PointerWrapper<double>, const py::object&)>;

   using PythonUserModel = UserModel<std::optional<Objective>, std::optional<ObjectiveGradient>, std::optional<Constraints>,
      std::optional<Jacobian>, std::optional<JacobianOperator>, std::optional<JacobianTransposedOperator>,
      std::optional<Hessian>, std::optional<HessianOperator>, std::optional<std::vector<double>>, std::optional<py::object>>;
} // namespace

#endif // UNO_UNOPY_H