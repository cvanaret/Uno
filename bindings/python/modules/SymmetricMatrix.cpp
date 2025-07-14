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
         .def("insert", [](PythonSymmetricMatrix& self, size_t row_index, size_t column_index, double element) {
            self->insert(row_index, column_index, element);
         }, py::arg("row_index"), py::arg("column_index"), py::arg("element"), "Insert an element");
   }
} // namespace