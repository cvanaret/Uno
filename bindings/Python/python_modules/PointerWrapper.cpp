// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "tools/PointerWrapper.hpp"

namespace py = pybind11;

namespace uno {
   void define_PointerWrapper(py::module& module) {
      py::class_<PointerWrapper<double>>(module, "PointerToDouble")
         // accessor/setter
         .def("__getitem__", [](const PointerWrapper<double> self, size_t index) {
            return self[index];
         })
         .def("__setitem__", [](PointerWrapper<double> self, size_t index, double value) {
            self[index] = value;
         });

      py::class_<PointerWrapper<const double>>(module, "ConstPointerToDouble")
         // accessor/setter
         .def("__getitem__", [](const PointerWrapper<const double> self, size_t index) {
            return self[index];
         });

      py::class_<PointerWrapper<void>>(module, "PointerToVoid");
   }
} // namespace