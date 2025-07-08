// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "linear_algebra/SymmetricMatrix.hpp"

namespace py = pybind11;

namespace uno {
   void define_SymmetricMatrix(py::module& module) {
      // SymmetricMatrix class (instantiated with size_t indices and double elements)
      py::class_<SymmetricMatrix<size_t, double>>(module, "SymmetricMatrix")
         // abstract class, no constructor
         // methods
         .def("insert", &SymmetricMatrix<size_t, double>::insert, py::arg("element"), py::arg("row_index"),
            py::arg("column_index"), "Insert an element")
         .def("finalize_column", &SymmetricMatrix<size_t, double>::finalize_column, py::arg("column_index"),
            "Finalize a given column");
   }
} // namespace