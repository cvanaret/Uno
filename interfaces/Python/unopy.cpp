// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "unopy.hpp"
#include "Uno.hpp"

namespace py = pybind11;

namespace uno {
   void define_Model(py::module& module);
   void define_Result(py::module& module);
   void define_UnoSolver(py::module& module);
   void define_Constants(py::module& module);

   // unopy module definition
   PYBIND11_MODULE(unopy, module) {
      std::string description = "Python bindings to Uno ";
      description.append(Uno::current_version());
      description.append(", a solver for nonlinearly constrained optimization");
      module.doc() = description;

      module.def("current_uno_version", &Uno::current_version);

      // raised by a callback (model evaluation or solver callback) to stop Uno immediately with USER_TERMINATION
      static py::exception<UserTermination> user_termination(module, "UserTermination", PyExc_Exception);
      user_termination.attr("__doc__") = "Raise it in a callback (model evaluation or solver callback) to stop Uno "
         "immediately. optimize() then returns the last accepted iterate with OptimizationStatus.USER_TERMINATION.";
      user_termination_type() = user_termination;

      define_Model(module);
      define_Result(module);
      define_UnoSolver(module);
      define_Constants(module);
   }
} // namespace