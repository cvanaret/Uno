// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include "../classes/UnoSolverWrapper.hpp"
#include "options/Options.hpp"

namespace py = pybind11;

namespace uno {
   void define_UnoSolver(py::module& module) {
      py::class_<UnoSolverWrapper>(module, "UnoSolver")
         // constructor
         .def(py::init<bool, const Options&>(), py::arg("constrained_model"), py::arg("options"), "Constructor")
         // methods
         .def("solve", &UnoSolverWrapper::solve,
            /*
            py::arg("number_variables"), py::arg("number_constraints"), py::arg("evaluate_objective"), py::arg("evaluate_constraints"),
            py::arg("evaluate_objective_gradient"), py::arg("evaluate_jacobian"), py::arg("evaluate_lagrangian_hessian"),
            py::arg("options"),
            */
            "Solve an optimization model with the Uno solver");
   }
} // namespace