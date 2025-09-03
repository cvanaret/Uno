// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <string>
#include "../unopy.hpp"

namespace py = pybind11;

namespace uno {
   void define_Model(py::module& module) {
      py::class_<PythonUserModel>(module, "Model")
         // constructor
         .def(py::init<char, int32_t, double*, double*, int32_t>(), "Constructor");
      /*
      [](char problem_type, int32_t number_variables, std::vector<double>& variables_lower_bounds,
      std::vector<double>& variables_upper_bounds, int32_t base_indexing) {
      return PythonUserModel(problem_type, number_variables, variables_lower_bounds.data(),
      variables_upper_bounds.data(), base_indexing);
      }
      */
   }
} // namespace