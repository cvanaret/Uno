// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include "optimization/Result.hpp"

namespace py = pybind11;

namespace uno {
   void define_Result(py::module& module) {
      py::class_<Result>(module, "Result")
      // attributes
      .def_readonly("optimization_status", &Result::optimization_status)
      .def_readonly("solution_status", &Result::solution_status)
      .def_readonly("solution_objective", &Result::solution_objective)
      .def_readonly("solution_primal_feasibility", &Result::solution_primal_feasibility)
      .def_readonly("solution_dual_feasibility", &Result::solution_dual_feasibility)
      .def_readonly("solution_complementarity", &Result::solution_complementarity)
      .def_readonly("primal_solution", &Result::primal_solution)
      .def_readonly("constraint_dual_solution", &Result::constraint_dual_solution)
      .def_readonly("lower_bound_dual_solution", &Result::lower_bound_dual_solution)
      .def_readonly("upper_bound_dual_solution", &Result::upper_bound_dual_solution)
      .def_readonly("number_iterations", &Result::number_iterations)
      .def_readonly("cpu_time", &Result::cpu_time)
      .def_readonly("number_objective_evaluations", &Result::number_objective_evaluations)
      .def_readonly("number_constraint_evaluations", &Result::number_constraint_evaluations)
      .def_readonly("number_objective_gradient_evaluations", &Result::number_objective_gradient_evaluations)
      .def_readonly("number_jacobian_evaluations", &Result::number_jacobian_evaluations)
      .def_readonly("number_hessian_evaluations", &Result::number_hessian_evaluations)
      .def_readonly("number_subproblems_solved", &Result::number_subproblems_solved);

      py::enum_<OptimizationStatus>(module, "OptimizationStatus")
      .value("SUCCESS", OptimizationStatus::SUCCESS)
      .value("ITERATION_LIMIT", OptimizationStatus::ITERATION_LIMIT)
      .value("TIME_LIMIT", OptimizationStatus::TIME_LIMIT)
      .value("EVALUATION_ERROR", OptimizationStatus::EVALUATION_ERROR)
      .value("ALGORITHMIC_ERROR", OptimizationStatus::ALGORITHMIC_ERROR);

      py::enum_<SolutionStatus>(module, "SolutionStatus")
      .value("NOT_OPTIMAL", SolutionStatus::NOT_OPTIMAL)
      .value("FEASIBLE_KKT_POINT", SolutionStatus::FEASIBLE_KKT_POINT)
      .value("FEASIBLE_FJ_POINT", SolutionStatus::FEASIBLE_FJ_POINT)
      .value("INFEASIBLE_STATIONARY_POINT", SolutionStatus::INFEASIBLE_STATIONARY_POINT)
      .value("FEASIBLE_SMALL_STEP", SolutionStatus::FEASIBLE_SMALL_STEP)
      .value("INFEASIBLE_SMALL_STEP", SolutionStatus::INFEASIBLE_SMALL_STEP)
      .value("UNBOUNDED", SolutionStatus::UNBOUNDED);
   }
} // namespace