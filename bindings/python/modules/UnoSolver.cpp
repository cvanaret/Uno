// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "Uno.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"

namespace py = pybind11;

namespace uno {
   void define_UnoSolver(py::module& module) {
      // Uno class
      py::class_<Uno>(module, "UnoSolver")
         // constructor
         .def(py::init<size_t, const Options&>(), py::arg("number_constraints"), py::arg("options"), "Constructor")
         // methods
         .def("solve", py::overload_cast<const Model&, Iterate&, const Options&>(&Uno::solve),
            py::arg("model"), py::arg("iterate"), py::arg("options"),
            "Solve an optimization model with the Uno solver");
      ;
   }
} // namespace