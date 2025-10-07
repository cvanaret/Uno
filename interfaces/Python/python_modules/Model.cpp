// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include "../unopy.hpp"
#include "symbolic/Range.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(uno::Vector<double>)

namespace uno {
   void define_Model(py::module& module) {
      py::class_<PythonUserModel>(module, "Model")
      // constructor
      // https://stackoverflow.com/a/62310838/16037994
      .def(py::init<>([](char problem_type, int32_t number_variables, const std::vector<double>& variables_lower_bounds,
            const std::vector<double>& variables_upper_bounds, int32_t base_indexing) {
         PythonUserModel model(problem_type, number_variables, base_indexing);
         // copy the bounds internally
         const size_t unsigned_number_variables = static_cast<size_t>(number_variables);
         model.variables_lower_bounds = std::vector<double>(unsigned_number_variables);
         model.variables_upper_bounds = std::vector<double>(unsigned_number_variables);
         for (size_t variable_index: Range(unsigned_number_variables)) {
            (*model.variables_lower_bounds)[variable_index] = variables_lower_bounds[variable_index];
            (*model.variables_upper_bounds)[variable_index] = variables_upper_bounds[variable_index];
         }
         return model;
      }), "Constructor")

      // methods
      .def("set_objective", [](PythonUserModel& user_model, int32_t optimization_sense, const Objective& objective_function,
            const ObjectiveGradient& objective_gradient) {
         if (optimization_sense != UNO_MINIMIZE && optimization_sense != UNO_MAXIMIZE) {
            throw std::runtime_error("Please specify a valid objective sense.");
         }
         user_model.optimization_sense = optimization_sense;
         user_model.objective_function = objective_function;
         user_model.objective_gradient = objective_gradient;
      })

      .def("set_constraints", [](PythonUserModel& user_model, int32_t number_constraints, const Constraints& constraint_functions,
            std::vector<double>& constraints_lower_bounds, std::vector<double>& constraints_upper_bounds, int32_t number_jacobian_nonzeros,
            const std::vector<int32_t>& jacobian_row_indices, const std::vector<int32_t>& jacobian_column_indices,
            const Jacobian& constraint_jacobian) {
         if (number_constraints <= 0) {
            throw std::runtime_error("Please specify a positive number of constraints.");
         }
         user_model.number_constraints = number_constraints;
         user_model.constraint_functions = constraint_functions;
         user_model.constraints_lower_bounds = constraints_lower_bounds;
         user_model.constraints_upper_bounds = constraints_upper_bounds;
         user_model.number_jacobian_nonzeros = number_jacobian_nonzeros;
         // copy the Jacobian sparsity to allow the calling code to dispose of its vectors
         user_model.jacobian_row_indices.resize(static_cast<size_t>(number_jacobian_nonzeros));
         user_model.jacobian_column_indices.resize(static_cast<size_t>(number_jacobian_nonzeros));
         for (size_t index: Range(static_cast<size_t>(number_jacobian_nonzeros))) {
            user_model.jacobian_row_indices[index] = jacobian_row_indices[index];
            user_model.jacobian_column_indices[index] = jacobian_column_indices[index];
         }
         user_model.constraint_jacobian = constraint_jacobian;
      })

      .def("set_lagrangian_hessian", [](PythonUserModel& user_model, int32_t number_hessian_nonzeros, char hessian_triangular_part,
            const std::vector<int32_t>& hessian_row_indices, const std::vector<int32_t>& hessian_column_indices,
            const Hessian& lagrangian_hessian, double lagrangian_sign_convention) {
         if (number_hessian_nonzeros <= 0) {
            throw std::runtime_error("Please specify a positive number of Lagrangian Hessian nonzeros.");
         }
         if (hessian_triangular_part != UNO_LOWER_TRIANGLE && hessian_triangular_part != UNO_UPPER_TRIANGLE) {
            throw std::runtime_error("Please specify a correct Hessian triangle.");
         }
         if (lagrangian_sign_convention != UNO_MULTIPLIER_NEGATIVE && lagrangian_sign_convention != UNO_MULTIPLIER_POSITIVE) {
            throw std::runtime_error("Please specify a correct Lagrangian sign convention.");
         }
         
         user_model.number_hessian_nonzeros = number_hessian_nonzeros;
         // copy the Hessian sparsity to allow the calling code to dispose of its vectors
         // from now on, we only maintain the lower triangle
         user_model.hessian_row_indices.resize(static_cast<size_t>(number_hessian_nonzeros));
         user_model.hessian_column_indices.resize(static_cast<size_t>(number_hessian_nonzeros));
         const bool lower_triangle = (hessian_triangular_part == UNO_LOWER_TRIANGLE);
         for (size_t index: Range(static_cast<size_t>(number_hessian_nonzeros))) {
            user_model.hessian_row_indices[index] = lower_triangle ? hessian_row_indices[index] : hessian_column_indices[index];
            user_model.hessian_column_indices[index] = lower_triangle ? hessian_column_indices[index] : hessian_row_indices[index];
         }
         user_model.hessian_triangular_part = UNO_LOWER_TRIANGLE;
         user_model.lagrangian_hessian = lagrangian_hessian;
         user_model.lagrangian_sign_convention = lagrangian_sign_convention;
      })

      .def("set_jacobian_operator", [](PythonUserModel& user_model, const JacobianOperator& jacobian_operator) {
         user_model.jacobian_operator = jacobian_operator;
      })

      .def("set_jacobian_transposed_operator", [](PythonUserModel& user_model, const JacobianTransposedOperator& jacobian_transposed_operator) {
         user_model.jacobian_transposed_operator = jacobian_transposed_operator;
      })

      .def("set_lagrangian_hessian_operator", [](PythonUserModel& user_model, int32_t number_hessian_nonzeros,
            const HessianOperator& lagrangian_hessian_operator, double lagrangian_sign_convention) {
         if (number_hessian_nonzeros <= 0) {
            throw std::runtime_error("Please specify a positive number of Lagrangian Hessian nonzeros.");
         }
         if (lagrangian_sign_convention != UNO_MULTIPLIER_NEGATIVE && lagrangian_sign_convention != UNO_MULTIPLIER_POSITIVE) {
            throw std::runtime_error("Please specify a correct Lagrangian sign convention.");
         }
         
         user_model.number_hessian_nonzeros = number_hessian_nonzeros;
         user_model.lagrangian_hessian_operator = lagrangian_hessian_operator;
         user_model.lagrangian_sign_convention = lagrangian_sign_convention;
      })

      .def("set_initial_primal_iterate", [](PythonUserModel& user_model, std::vector<double>& initial_primal_iterate) {
         user_model.initial_primal_iterate = initial_primal_iterate;
      })

      .def("set_initial_dual_iterate", [](PythonUserModel& user_model, std::vector<double>& initial_dual_iterate) {
         user_model.initial_dual_iterate = initial_dual_iterate;
      })

      .def("set_user_data", [](PythonUserModel& user_model, const py::object& user_data) {
         user_model.user_data = user_data;
      });
   }
} // namespace