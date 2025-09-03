// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <string>
#include "../cpp_classes/UnoSolver.hpp"

namespace py = pybind11;

namespace uno {
   void define_UnoSolver(py::module& module) {
      py::class_<UnoSolver>(module, "UnoSolver")
         // constructor
      .def(py::init<>(), "Constructor")
         // methods
         .def("set_option", [](UnoSolver& solver, const std::string& option_name, const std::string& option_value) {
            solver.options.set(option_name, option_value);
         }, py::arg("option_name"), py::arg("option_value"))
         .def("optimize", &UnoSolver::optimize, py::arg("model"), "Optimize an optimization model with the Uno solver");
   }
} // namespace