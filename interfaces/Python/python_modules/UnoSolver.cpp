// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <memory>
#include <string>
#include "../cpp_classes/PythonModel.hpp"
#include "../cpp_classes/PythonStreamBuffer.hpp"
#include "../cpp_classes/UnoSolverWrapper.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"

namespace py = pybind11;

namespace uno {
   void define_UnoSolver(py::module& module) {
      py::class_<UnoSolverWrapper>(module, "UnoSolver")
      // constructor
      .def(py::init<>(), "Constructor")

      // methods
      .def("set_option", [](UnoSolverWrapper& solver, const std::string& option_name, bool option_value) {
         const auto type = Options::option_types.find(option_name);
         if (type != Options::option_types.end()) { // key found
            if (type->second == OptionType::BOOL) { // correct type
               solver.options.set_bool(option_name, option_value);
            }
            else { // incorrect type
               throw py::type_error(option_name + " is not of type bool");
            }
         }
         else { // key not found
            throw py::key_error(option_name + " does not exist");
         }
      }, py::arg("option_name"), py::arg("option_value"))

      .def("set_option", [](UnoSolverWrapper& solver, const std::string& option_name, uno_int option_value) {
         const auto type = Options::option_types.find(option_name);
         if (type != Options::option_types.end()) { // key found
            if (type->second == OptionType::INTEGER) { // correct type
               solver.options.set_integer(option_name, option_value);
            }
            else { // incorrect type
               throw py::type_error(option_name + " is not of type int");
            }
         }
         else { // key not found
            throw py::key_error(option_name + " does not exist");
         }
      }, py::arg("option_name"), py::arg("option_value"))

      .def("set_option", [](UnoSolverWrapper& solver, const std::string& option_name, double option_value) {
         const auto type = Options::option_types.find(option_name);
         if (type != Options::option_types.end()) { // key found
            if (type->second == OptionType::DOUBLE) { // correct type
               solver.options.set_double(option_name, option_value);
            }
            else { // incorrect type
               throw py::type_error(option_name + " is not of type double");
            }
         }
         else { // key not found
            throw py::key_error(option_name + " does not exist");
         }
      }, py::arg("option_name"), py::arg("option_value"))

      .def("set_option", [](UnoSolverWrapper& solver, const std::string& option_name, const std::string& option_value) {
         const auto type = Options::option_types.find(option_name);
         if (type != Options::option_types.end()) { // key found
            if (type->second == OptionType::STRING) { // correct type
               solver.options.set_string(option_name, option_value);
            }
            else { // incorrect type
               throw py::type_error(option_name + " is not of type string");
            }
         }
         else { // key not found
            throw py::key_error(option_name + " does not exist");
         }
      }, py::arg("option_name"), py::arg("option_value"))

      .def("set_logger_stream", &UnoSolverWrapper::set_logger_stream)

      .def("set_preset", [](UnoSolverWrapper& solver, const std::string& preset_name) {
         Presets::set(solver.options, preset_name);
      }, py::arg("preset_name"))

      .def("optimize", [](UnoSolverWrapper& solver, const PythonUserModel& user_model) {
         return solver.optimize(user_model);
      }, py::arg("model"), "Optimize an optimization model with the Uno solver");
   }
} // namespace