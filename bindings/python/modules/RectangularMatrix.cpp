// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "tools/PointerWrapper.hpp"

namespace py = pybind11;

namespace uno {
   using PythonRectangularMatrix = PointerWrapper<RectangularMatrix<double>>;
   using PythonSparseVector = PointerWrapper<SparseVector<double>>;

   void define_RectangularMatrix(py::module& module) {
      // RectangularMatrix class (instantiated with double elements)
      py::class_<PythonRectangularMatrix>(module, "RectangularMatrix")
         .def("__getitem__", [](PythonRectangularMatrix& self, size_t row_index) -> PythonSparseVector {
            return wrap_pointer(&(*self)[row_index]);
         }, py::arg("row_index"), "Returns a row", py::return_value_policy::reference_internal);
   }
} // namespace