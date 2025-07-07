// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <string>
#include "model/Model.hpp"

namespace py = pybind11;

namespace uno {
   void define_Model(py::module& module) {
      // Model class
      py::class_<Model>(module, "Model")
         // constructor
         .def(py::init<std::string, size_t, size_t, double>(), py::arg("name"), py::arg("number_variables"),
            py::arg("number_constraints"), py::arg("objective_sign"), "Constructor")
      ;
   }
} // namespace