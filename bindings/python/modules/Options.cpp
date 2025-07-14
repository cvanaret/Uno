// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"

namespace py = pybind11;

namespace uno {
   void define_Options(py::module& module) {
      // Options class
      py::class_<Options>(module, "Options")
         // constructor
         .def(py::init<bool>(), py::arg("are_default_options"), "Constructor")
         // accessor/setter
         .def("__getitem__", &Options::get_string, py::arg("key"), "Read an option")
         .def("__setitem__", [](Options& options, const std::string& key, const std::string& value) {
            options[key] = value;
         }, py::arg("key"), py::arg("value"), "Set an option")
         // string representation
         .def("__repr__", [](const Options& options) {
            return options.to_string();
         });
      // get the default options
      module.def("get_default_options", []() {
         Options options = DefaultOptions::load();
         // determine the default solvers based on the available libraries
         Options solvers_options = DefaultOptions::determine_solvers();
         options.overwrite_with(solvers_options);
         return options;
      }, "create default options");
      // set a preset
      module.def("set_preset", [](Options& options, const std::string& preset_name) {
         Presets::set(options, preset_name);
      }, "set a preset");
   }
} // namespace