// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "unopy.hpp"
#include "Uno.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(uno::Vector<double>);

namespace uno {
   void define_PythonUserModel(py::module& module);
   void define_UnoSolver(py::module& module);

   // unopy module definition
   PYBIND11_MODULE(unopy, module) {
      std::string description = "Python bindings to Uno ";
      description.append(Uno::current_version());
      description.append(", a solver for nonlinearly constrained optimization");
      module.doc() = description;

      define_PythonUserModel(module);
      define_UnoSolver(module);
   }
} // namespace