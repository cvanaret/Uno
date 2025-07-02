#include <pybind11/pybind11.h>
#include <string>
#include "model/Model.hpp"

namespace py = pybind11;

namespace uno {
   void define_Model(py::module& module) {
      // Model class
      py::class_<Model>(module, "Model")
         // constructor
         .def(py::init<std::string, size_t, size_t, double>(), py::arg("std::string name"), py::arg("size_t number_variables"),
            py::arg("size_t number_constraints"), py::arg("double objective_sign"), "Constructor")
      ;
   }
} // namespace
