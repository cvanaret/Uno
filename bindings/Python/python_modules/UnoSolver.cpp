// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include "../cpp_classes/PythonModel.hpp"
#include "../cpp_classes/UnoSolverWrapper.hpp"
#include "options/Presets.hpp"

namespace py = pybind11;

namespace uno {
   void define_UnoSolver(py::module& module) {
      py::class_<UnoSolverWrapper>(module, "UnoSolver")
         // constructor
         .def(py::init<>(), "Constructor")

         // methods
         .def("set_option", [](UnoSolverWrapper& solver, const std::string& option_name, const std::string& option_value) {
            solver.options.set(option_name, option_value);
         }, py::arg("option_name"), py::arg("option_value"))

         .def("set_preset", [](UnoSolverWrapper& solver, const std::string& preset_name) {
            Presets::set(solver.options, preset_name);
         }, py::arg("preset_name"))

         .def("optimize", [](UnoSolverWrapper& solver, const PythonUserModel& user_model) {
            return solver.optimize(user_model);
         }, py::arg("model"), "Optimize an optimization model with the Uno solver");
   }
} // namespace