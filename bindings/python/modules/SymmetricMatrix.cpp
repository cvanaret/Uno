// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "linear_algebra/SymmetricMatrix.hpp"

namespace py = pybind11;

namespace uno {
   void define_SymmetricMatrix(py::module& module) {
      // RectangularMatrix class (instantiated with double elements)
      py::class_<SymmetricMatrix<size_t, double>>(module, "SymmetricMatrix")
         // abstract class, no constructor
         // methods
         //.def("get_number_rows", &RectangularMatrix<double>::get_number_rows, "Number of rows")
         //.def("get_number_columns", &RectangularMatrix<double>::get_number_columns, "Number of columns")
         .def("insert", &SymmetricMatrix<size_t, double>::insert, py::arg("element"), py::arg("row_index"),
            py::arg("column_index"), "Insert an element");
         /*
         .def("__getitem__", [](RectangularMatrix<double>& matrix, size_t row_index) -> SparseVector<double>& {
            return matrix[row_index];
         }, py::arg("row_index"), "Returns a row", py::return_value_policy::reference_internal);
         */
   }
} // namespace