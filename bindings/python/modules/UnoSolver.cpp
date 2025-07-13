// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <functional>
#include "Uno.hpp"
#include "options/Options.hpp"

namespace py = pybind11;

namespace uno {
   void define_UnoSolver(py::module& module) {
      py::class_<Uno>(module, "UnoSolver")
         // constructor
         .def(py::init<bool, const Options&>(), py::arg("constrained_model"), py::arg("options"), "Constructor")
         // methods
         .def("solve", py::overload_cast<size_t, size_t, const Uno::objective_function_type&,
               const Uno::constraint_functions_type&, const Options&>(&Uno::solve),
            py::arg("number_variables"), py::arg("number_constraints"), py::arg("objective"), py::arg("constraints"),
            py::arg("options"), "Solve an optimization model with the Uno solver");
      ;
   }
} // namespace