// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "linear_algebra/SymmetricMatrix.hpp"
#include "tools/PointerWrapper.hpp"

namespace py = pybind11;

namespace uno {
   using PythonSymmetricMatrix = PointerWrapper<SymmetricMatrix<size_t, double>>;

   void define_SymmetricMatrix(py::module& module) {
      // SymmetricMatrix class (instantiated with size_t indices and double elements)
      py::class_<PythonSymmetricMatrix>(module, "SymmetricMatrix")
         // methods
         .def("insert", [](PythonSymmetricMatrix& self, double element, size_t row_index, size_t column_index) {
            self->insert(element, row_index, column_index);
         }, py::arg("element"), py::arg("row_index"), py::arg("column_index"), "Insert an element");
   }
} // namespace