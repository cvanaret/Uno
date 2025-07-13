// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "linear_algebra/Vector.hpp"
#include "tools/PointerWrapper.hpp"

namespace py = pybind11;

namespace uno {
   using PythonVector = PointerWrapper<Vector<double>>;

   void define_Vector(py::module& module) {
      // Vector class (instantiated with double elements)
      py::class_<PythonVector>(module, "Vector")
         // constructor
         //.def(py::init<size_t>(), py::arg("capacity"), "Constructor")
         // methods
         .def("size", [](const PythonVector self) {
            return self->size();
         }, "Number of elements of the vector")
         // iterator
         .def("__iter__", [](const PythonVector self) {
            return py::make_iterator(self->begin(), self->end());
         }, "Iterator on the (index, value) pairs", py::keep_alive<0, 1>())
         // accessor/setter
         .def("__getitem__", [](const PythonVector self, size_t index) {
            assert(self.operator->() != nullptr && "pybind11's Vector __getitem__: dereferencing NULL pointer");
            return (*self)[index];
         })
         .def("__setitem__", [](PythonVector self, size_t index, double value) {
            assert(self.operator->() != nullptr && "pybind11's Vector __setitem__: dereferencing NULL pointer");
            (*self)[index] = value;
         });
   }
} // namespace