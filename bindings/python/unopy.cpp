// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "linear_algebra/Vector.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(uno::Vector<double>);

namespace uno {
   void define_SparseVector(py::module& module);
   void define_RectangularMatrix(py::module& module);
   void define_SymmetricMatrix(py::module& module);
   void define_Options(py::module& module);
   void define_UnoSolver(py::module& module);
   void define_Vector(py::module& module);

   // unopy module definition
   PYBIND11_MODULE(unopy, module) {
      module.doc() = "Python bindings to the solver Uno for nonlinearly constrained optimization";
      
      define_SparseVector(module);
      define_RectangularMatrix(module);
      define_SymmetricMatrix(module);
      define_Options(module);
      define_UnoSolver(module);
      define_Vector(module);
   }
} // namespace