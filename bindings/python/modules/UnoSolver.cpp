// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <functional>
#include <vector>
#include "Uno.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"

namespace py = pybind11;

using objective_type = const std::function<double(const std::vector<double>&)>&;

namespace uno {
   void define_UnoSolver(py::module& module) {
      // Uno class
      py::class_<Uno>(module, "UnoSolver")
         // constructor
         .def(py::init<bool, const Options&>(), py::arg("constrained_model"), py::arg("options"), "Constructor")
         // methods
         .def("solve", py::overload_cast<size_t, size_t, objective_type>(&Uno::solve),
            py::arg("number_variables"), py::arg("number_constraints"), py::arg("objective"),"Dummy call")
         .def("solve", py::overload_cast<const Model&, Iterate&, const Options&>(&Uno::solve),
            py::arg("model"), py::arg("iterate"), py::arg("options"),
            "Solve an optimization model with the Uno solver");
      ;
   }
} // namespace