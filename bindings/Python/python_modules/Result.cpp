// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "optimization/Result.hpp"

namespace py = pybind11;

namespace uno {
   void define_Result(py::module& module) {
      py::class_<Result>(module, "Result")
      // attributes
      //.def_readonly("primal_solution", &Result::solution.primals)
      .def_readonly("optimization_status", &Result::optimization_status)
      .def_readonly("number_iterations", &Result::number_iterations)
      .def_readonly("cpu_time", &Result::cpu_time)
      .def_readonly("objective_evaluations", &Result::number_objective_evaluations)
      .def_readonly("constraint_evaluations", &Result::number_constraint_evaluations)
      .def_readonly("objective_gradient_evaluations", &Result::number_objective_gradient_evaluations)
      .def_readonly("jacobian_evaluations", &Result::number_jacobian_evaluations)
      .def_readonly("hessian_evaluations", &Result::number_hessian_evaluations)
      .def_readonly("number_subproblems_solved", &Result::number_subproblems_solved);
   }
} // namespace