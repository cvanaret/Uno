// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <utility>
#include <stdexcept>
#include "../unopy.hpp"
#include "tools/Infinity.hpp"

namespace py = pybind11;

namespace uno {
   void define_Model(py::module& module) {
      py::class_<PythonUserModel>(module, "Model")
      // constructor with bounds
      // https://stackoverflow.com/a/62310838/16037994
      .def(py::init<>([](const std::string& problem_type, uno_int number_variables, ArrayLike<double> variables_lower_bounds,
            ArrayLike<double> variables_upper_bounds, uno_int base_indexing) {
         PythonUserModel model(problem_type.data(), number_variables, base_indexing);
         const size_t unsigned_number_variables = static_cast<size_t>(number_variables);
         model.variables_lower_bounds = to_vector(variables_lower_bounds, unsigned_number_variables, "variables_lower_bounds");
         model.variables_upper_bounds = to_vector(variables_upper_bounds, unsigned_number_variables, "variables_upper_bounds");
         model.initial_primal_iterate.resize(static_cast<size_t>(number_variables), 0.);
         return model;
      }), "Constructor")

      // constructor without bounds
      .def(py::init<>([](const std::string& problem_type, uno_int number_variables, uno_int base_indexing) {
         PythonUserModel model(problem_type.data(), number_variables, base_indexing);
         const size_t unsigned_number_variables = static_cast<size_t>(number_variables);
         model.variables_lower_bounds.resize(unsigned_number_variables, -Inf);
         model.variables_upper_bounds.resize(unsigned_number_variables, Inf);
         model.initial_primal_iterate.resize(unsigned_number_variables, 0.);
         return model;
      }), "Constructor")

      // methods
      .def("set_name", [](PythonUserModel& user_model, std::string name) {
         user_model.name = std::move(name);
      })

      .def("set_variables_lower_bounds", [](PythonUserModel& user_model, ArrayLike<double> variables_lower_bounds) {
         user_model.variables_lower_bounds = to_vector(variables_lower_bounds, static_cast<size_t>(user_model.number_variables), "variables_lower_bounds");
      })

      .def("set_variables_upper_bounds", [](PythonUserModel& user_model, ArrayLike<double> variables_upper_bounds) {
         user_model.variables_upper_bounds = to_vector(variables_upper_bounds, static_cast<size_t>(user_model.number_variables), "variables_upper_bounds");
      })

      .def("set_variable_lower_bound", [](PythonUserModel& user_model, uno_int variable_index, double lower_bound) {
         if (variable_index < 0 || variable_index >= user_model.number_variables) {
            throw std::invalid_argument("Please specify a valid index.");
         }
         user_model.variables_lower_bounds[static_cast<size_t>(variable_index)] = lower_bound;
      })

      .def("set_variable_upper_bound", [](PythonUserModel& user_model, uno_int variable_index, double upper_bound) {
         if (variable_index < 0 || variable_index >= user_model.number_variables) {
            throw std::invalid_argument("Please specify a valid index.");
         }
         user_model.variables_upper_bounds[static_cast<size_t>(variable_index)] = upper_bound;
      })

      .def("set_objective", [](PythonUserModel& user_model, uno_int optimization_sense, Objective objective_function,
            ObjectiveGradient objective_gradient) {
         if (optimization_sense != UNO_MINIMIZE && optimization_sense != UNO_MAXIMIZE) {
            throw std::invalid_argument("Please specify a valid objective sense.");
         }
         user_model.optimization_sense = optimization_sense;
         user_model.objective_function = std::move(objective_function);
         user_model.objective_gradient = std::move(objective_gradient);
      })

      .def("set_constraints", [](PythonUserModel& user_model, uno_int number_constraints, Constraints constraint_functions,
            ArrayLike<double> constraints_lower_bounds, ArrayLike<double> constraints_upper_bounds,
            uno_int number_jacobian_nonzeros, ArrayLike<uno_int> jacobian_row_indices,
            ArrayLike<uno_int> jacobian_column_indices, Jacobian jacobian) {
         if (number_constraints < 0) {
            throw std::invalid_argument("The number of constraints cannot be negative.");
         }
         if (number_jacobian_nonzeros < 0) {
            throw std::invalid_argument("The number of Jacobian nonzeros cannot be negative.");
         }
         const size_t m = static_cast<size_t>(number_constraints), nnz = static_cast<size_t>(number_jacobian_nonzeros);
         // convert everything before modifying the model, so that it is left untouched if an argument is invalid
         auto lower_bounds = to_vector(constraints_lower_bounds, m, "constraints_lower_bounds");
         auto upper_bounds = to_vector(constraints_upper_bounds, m, "constraints_upper_bounds");
         auto row_indices = to_vector(jacobian_row_indices, nnz, "jacobian_row_indices");
         auto column_indices = to_vector(jacobian_column_indices, nnz, "jacobian_column_indices");
         user_model.number_constraints = number_constraints;
         user_model.constraint_functions = std::move(constraint_functions);
         user_model.constraints_lower_bounds = std::move(lower_bounds);
         user_model.constraints_upper_bounds = std::move(upper_bounds);
         user_model.number_jacobian_nonzeros = number_jacobian_nonzeros;
         user_model.jacobian_row_indices = std::move(row_indices);
         user_model.jacobian_column_indices = std::move(column_indices);
         user_model.jacobian = std::move(jacobian);
         user_model.initial_dual_iterate.resize(static_cast<size_t>(number_constraints), 0.);
      })

      .def("set_constraints_lower_bounds", [](PythonUserModel& user_model, ArrayLike<double> constraints_lower_bounds) {
         user_model.constraints_lower_bounds = to_vector(constraints_lower_bounds, static_cast<size_t>(user_model.number_constraints), "constraints_lower_bounds");
      })

      .def("set_constraints_upper_bounds", [](PythonUserModel& user_model, ArrayLike<double> constraints_upper_bounds) {
         user_model.constraints_upper_bounds = to_vector(constraints_upper_bounds, static_cast<size_t>(user_model.number_constraints), "constraints_upper_bounds");
      })

      .def("set_constraint_lower_bound", [](PythonUserModel& user_model, uno_int constraint_index, double lower_bound) {
         if (constraint_index < 0 || constraint_index >= user_model.number_constraints) {
            throw std::invalid_argument("Please specify a valid index.");
         }
         user_model.constraints_lower_bounds[static_cast<size_t>(constraint_index)] = lower_bound;
      })

      .def("set_constraint_upper_bound", [](PythonUserModel& user_model, uno_int constraint_index, double upper_bound) {
         if (constraint_index < 0 || constraint_index >= user_model.number_constraints) {
            throw std::invalid_argument("Please specify a valid index.");
         }
         user_model.constraints_upper_bounds[static_cast<size_t>(constraint_index)] = upper_bound;
      })

      .def("set_lagrangian_hessian", [](PythonUserModel& user_model, uno_int number_hessian_nonzeros, char hessian_triangular_part,
            ArrayLike<uno_int> hessian_row_indices, ArrayLike<uno_int> hessian_column_indices, Hessian lagrangian_hessian) {
         if (number_hessian_nonzeros < 0) {
            throw std::invalid_argument("The number of Lagrangian Hessian nonzeros cannot be negative.");
         }
         if (hessian_triangular_part != UNO_LOWER_TRIANGLE && hessian_triangular_part != UNO_UPPER_TRIANGLE) {
            throw std::invalid_argument("Please specify a correct Hessian triangle.");
         }

         user_model.number_hessian_nonzeros = number_hessian_nonzeros;
         // if the Hessian is empty, the problem is an LP
         if (number_hessian_nonzeros == 0) {
            user_model.problem_type = ProblemType::LINEAR;
         }
         else {
            const size_t nnz = static_cast<size_t>(number_hessian_nonzeros);
            auto row_indices = to_vector(hessian_row_indices, nnz, "hessian_row_indices");
            auto column_indices = to_vector(hessian_column_indices, nnz, "hessian_column_indices");
            // we only maintain the lower triangle
            if (hessian_triangular_part == UNO_UPPER_TRIANGLE) {
               std::swap(row_indices, column_indices);
            }
            user_model.hessian_row_indices = std::move(row_indices);
            user_model.hessian_column_indices = std::move(column_indices);
            user_model.hessian_triangular_part = UNO_LOWER_TRIANGLE;
            user_model.lagrangian_hessian = std::move(lagrangian_hessian);
         }
      })

      .def("set_jacobian_operator", [](PythonUserModel& user_model, JacobianOperator jacobian_operator) {
         user_model.jacobian_operator = std::move(jacobian_operator);
      })

      .def("set_jacobian_transposed_operator", [](PythonUserModel& user_model, JacobianTransposedOperator jacobian_transposed_operator) {
         user_model.jacobian_transposed_operator = std::move(jacobian_transposed_operator);
      })

      .def("set_lagrangian_hessian_operator", [](PythonUserModel& user_model, HessianOperator lagrangian_hessian_operator) {
         user_model.lagrangian_hessian_operator = std::move(lagrangian_hessian_operator);
      })

      .def("set_lagrangian_sign_convention", [](PythonUserModel& user_model, uno_int lagrangian_sign_convention) {
         if (lagrangian_sign_convention != UNO_MULTIPLIER_NEGATIVE && lagrangian_sign_convention != UNO_MULTIPLIER_POSITIVE) {
            throw std::invalid_argument("Please specify a correct Lagrangian sign convention.");
         }
         user_model.lagrangian_sign_convention = lagrangian_sign_convention;
      })

      .def("set_initial_primal_iterate", [](PythonUserModel& user_model, ArrayLike<double> initial_primal_iterate) {
         user_model.initial_primal_iterate = to_vector(initial_primal_iterate, static_cast<size_t>(user_model.number_variables), "initial_primal_iterate");
      })

      .def("set_initial_dual_iterate", [](PythonUserModel& user_model, ArrayLike<double> initial_dual_iterate) {
         user_model.initial_dual_iterate = to_vector(initial_dual_iterate, static_cast<size_t>(user_model.number_constraints), "initial_dual_iterate");
      });
   }
} // namespace
