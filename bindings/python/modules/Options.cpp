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
         .def(py::init<>(), "Constructor")
         // accessor/setter
         .def("__getitem__", &Options::get_string, py::arg("key"), "Read an option")
         .def("__setitem__", [](Options& options, const std::string& key, const std::string& value) {
            options.set(key, value);
         }, py::arg("key"), py::arg("value"), "Set an option")
         // string representation
         .def("__repr__", [](const Options& /*options*/) {
            return ""; //options.to_string(); // TODO
         });
      // get the default options
      module.def("get_default_options", []() {
         Options options;
         DefaultOptions::load(options);
         return options;
      }, "create default options");
      // set a preset
      module.def("set_preset", [](Options& options, const std::string& preset_name) {
         const Options preset_options = Presets::get_preset_options(preset_name);
         options.overwrite_with(preset_options);
      }, "set a preset");
   }
} // namespace