// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include <pybind11/pybind11.h>
#include "linear_algebra/Vector.hpp"

namespace py = pybind11;

namespace uno {
   void define_Vector(py::module& module) {
      py::class_<Vector<double>>(module, "Vector")
         // methods
         .def("size", [](const Vector<double>& self) {
            return self.size();
         }, "Number of elements of the vector")
         // representation
         .def("__repr__", [](const Vector<double>& self) {
            std::string representation{};
            for (size_t index = 0; index < self.size(); ++index) {
               if (0 < index) {
                  representation.append(", ");
               }
               representation.append(std::to_string(self[index]));
            }
            return representation;
         })
         // iterator
         .def("__iter__", [](const Vector<double>& self) {
            return py::make_iterator(self.begin(), self.end());
         }, "Iterator on the (index, value) pairs", py::keep_alive<0, 1>())
         // accessor/setter
         .def("__getitem__", [](const Vector<double>& self, size_t index) {
            return self[index];
         }, py::return_value_policy::reference_internal)
         .def("__setitem__", [](Vector<double>& self, size_t index, double value) {
            self[index] = value;
         });
   }
} // namespace