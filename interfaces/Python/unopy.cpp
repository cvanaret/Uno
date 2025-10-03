// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "unopy.hpp"
#include "Uno.hpp"
#include "tools/PointerWrapper.hpp"

namespace py = pybind11;

//PYBIND11_MAKE_OPAQUE(uno::Vector<double>)

namespace uno {
   void define_Vector(py::module& module);
   void define_PointerWrapper(py::module& module);
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

      define_Vector(module);
      define_PointerWrapper(module);
      define_Model(module);
      define_Result(module);
      define_UnoSolver(module);
      define_Constants(module);
   }
} // namespace