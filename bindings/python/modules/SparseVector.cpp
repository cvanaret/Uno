// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "linear_algebra/SparseVector.hpp"
#include "tools/PointerWrapper.hpp"

namespace py = pybind11;

namespace uno {
   using PythonSparseVector = PointerWrapper<SparseVector<double>>;

   void define_SparseVector(py::module& module) {
      // SparseVector class (instantiated with double elements)
      py::class_<PythonSparseVector>(module, "SparseVector")
         //.def("size", &SparseVector<double>::size, "Number of elements")
         .def("insert", [](PythonSparseVector& self, size_t index, double value) {
            self->insert(index, value);
         }, py::arg("index"), py::arg("value"), "Insert an element");
         //.def("clear", &SparseVector<double>::clear, "Delete all elements")
         //.def("is_empty", &SparseVector<double>::is_empty, "True if the sparse vector is empty, false otherwise")
         // iterator
      /*
         .def("__iter__", [](const SparseVector<double>& vector) {
            return py::make_iterator(vector.begin(), vector.end());
         }, "Iterator on the (index, value) pairs", py::keep_alive<0, 1>())
         // string representation
         .def("__repr__", [](const SparseVector<double>& vector) {
            std::string representation = "Sparse vector:";
            for (const auto [index, value]: vector) {
               representation += "\n(" + std::to_string(index) + ", " + std::to_string(value) + ")";
            }
            return representation;
         });
         */
   }
} // namespace