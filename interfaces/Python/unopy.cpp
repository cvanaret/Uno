// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "unopy.hpp"
#include "Uno.hpp"
#include "ingredients/subproblem_solvers/Krylov/KrylovLoader.hpp"

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
      // fallback libkrylov, used when neither libkrylov_path nor UNO_KRYLOV_LIBRARY is set (set by __init__.py to the
      // library bundled in the wheel)
      module.def("_set_default_libkrylov_path", &set_default_krylov_library, py::arg("path"));
      module.def("libkrylov_error", &krylov_library_error,
         "Why libkrylov could not be loaded or is not functional (empty if no failure)");

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