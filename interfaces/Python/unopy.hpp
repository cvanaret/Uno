// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNOPY_H
#define UNO_UNOPY_H

#include <functional>
#include <optional>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "../UserModel.hpp"

namespace py = pybind11;

namespace uno {
   // alias for a typed, C-contiguous, non-owning numpy array view.
   // Use for INPUT: zero-copy if the source array already matches.
   template <typename T>
   using const_array = py::array_t<T, py::array::c_style | py::array::forcecast>;

   // for OUTPUT buffers owned by C++: wrap the raw pointer in a numpy array
   // with an explicit no-op capsule so Python knows not to free it.
   template <typename T>
   py::array_t<T> to_array(T* ptr, size_t n) {
      return py::array_t<T>(
          {static_cast<py::ssize_t>(n)}, // shape
          {sizeof(T)}, // strides
          ptr, // data (NOT owned)
          py::capsule(ptr, [](void*) {}) // no-op deleter: C++ owns this memory
      );
   }

   template <typename T>
   py::array_t<T> to_const_array(const T* ptr, size_t n) {
      auto array = py::array_t<double>(
          {static_cast<py::ssize_t>(n)},
          {sizeof(T)},
          const_cast<T*>(ptr),
          py::capsule(ptr, [](void*) {})
      );
      array.attr("flags").attr("writeable") = false;
      return array;
   }

   using Objective = std::function<int(const_array<double>, py::array_t<double>)>;

   using Constraints = std::function<int(const_array<double>, py::array_t<double>)>;

   using ObjectiveGradient = std::function<int(const_array<double>, py::array_t<double>)>;

   using Jacobian = std::function<int(const_array<double>, py::array_t<double>)>;

   using Hessian = std::function<int(const_array<double>, double, const_array<double>, py::array_t<double>)>;

   using JacobianOperator = std::function<int(const_array<double>, bool, const_array<double>, py::array_t<double>)>;

   using JacobianTransposedOperator = std::function<int(const_array<double>, bool, const_array<double>, py::array_t<double>)>;

   using HessianOperator = std::function<int(const_array<double>, bool, double, const_array<double>,
      const_array<double>, py::array_t<double>)>;

   using PythonUserModel = UserModel<std::optional<Objective>, std::optional<ObjectiveGradient>, std::optional<Constraints>,
      std::optional<Jacobian>, std::optional<JacobianOperator>, std::optional<JacobianTransposedOperator>,
      std::optional<Hessian>, std::optional<HessianOperator>, std::optional<std::vector<double>>, std::optional<py::object>>;
} // namespace

#endif // UNO_UNOPY_H