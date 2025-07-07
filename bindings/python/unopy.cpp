// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace uno {
   // individual classes
   void define_SparseVector(py::module& module);
   void define_RectangularMatrix(py::module& module);
   void define_SymmetricMatrix(py::module& module);
   void define_Options(py::module& module);
   //void define_Model(py::module& module);
   //void define_UnoSolver(py::module& module);

   // module definition
   PYBIND11_MODULE(unopy, module) {
      module.doc() = "Python bindings to the solver Uno for nonlinearly constrained optimization";
      
      define_SparseVector(module);
      define_RectangularMatrix(module);
      define_SymmetricMatrix(module);
      define_Options(module);
      //define_Model(module);
      //define_UnoSolver(module);
   }
} // namespace