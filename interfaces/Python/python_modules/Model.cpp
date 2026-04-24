// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include "../unopy.hpp"
#include "tools/Infinity.hpp"

namespace py = pybind11;

namespace uno {
   void define_Model(py::module& module) {
      py::class_<PythonUserModel>(module, "Model")
      // constructor with bounds
      // https://stackoverflow.com/a/62310838/16037994
      .def(py::init<>([](const std::string& problem_type, uno_int number_variables, std::vector<double> variables_lower_bounds,
            std::vector<double> variables_upper_bounds, uno_int base_indexing) {
         PythonUserModel model(problem_type.data(), number_variables, base_indexing);
         model.variables_lower_bounds = std::move(variables_lower_bounds);
         model.variables_upper_bounds = std::move(variables_upper_bounds);
         model.initial_primal_iterate.resize(static_cast<size_t>(number_variables), 0.);
         return model;
      }), "Constructor")

      // constructor without bounds
      .def(py::init<>([](const std::string& problem_type, uno_int number_variables, uno_int base_indexing) {
         PythonUserModel model(problem_type.data(), number_variables, base_indexing);
         const size_t unsigned_number_variables = static_cast<size_t>(number_variables);
         model.variables_lower_bounds.resize(unsigned_number_variables, -INF<double>);
         model.variables_upper_bounds.resize(unsigned_number_variables, INF<double>);
         model.initial_primal_iterate.resize(unsigned_number_variables, 0.);
         return model;
      }), "Constructor")

      // methods
      .def("set_variables_lower_bounds", [](PythonUserModel& user_model, std::vector<double> variables_lower_bounds) {
         user_model.variables_lower_bounds = std::move(variables_lower_bounds);
      })

      .def("set_variables_upper_bounds", [](PythonUserModel& user_model, std::vector<double> variables_upper_bounds) {
         user_model.variables_upper_bounds = std::move(variables_upper_bounds);
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
            std::vector<double> constraints_lower_bounds, std::vector<double> constraints_upper_bounds, uno_int number_jacobian_nonzeros,
            std::vector<uno_int> jacobian_row_indices, std::vector<uno_int> jacobian_column_indices, Jacobian jacobian) {
         if (number_constraints < 0) {
            throw std::invalid_argument("The number of constraints cannot be negative.");
         }
         if (constraints_lower_bounds.size() != static_cast<size_t>(number_constraints)) {
            throw std::invalid_argument("Dimension mismatch in constraints_lower_bounds.");
         }
         if (constraints_upper_bounds.size() != static_cast<size_t>(number_constraints)) {
            throw std::invalid_argument("Dimension mismatch in constraints_upper_bounds.");
         }
         if (jacobian_row_indices.size() != static_cast<size_t>(number_jacobian_nonzeros)) {
            throw std::invalid_argument("Dimension mismatch in jacobian_row_indices.");
         }
         if (jacobian_column_indices.size() != static_cast<size_t>(number_jacobian_nonzeros)) {
            throw std::invalid_argument("Dimension mismatch in jacobian_column_indices.");
         }
         user_model.number_constraints = number_constraints;
         user_model.constraint_functions = std::move(constraint_functions);
         user_model.constraints_lower_bounds = std::move(constraints_lower_bounds);
         user_model.constraints_upper_bounds = std::move(constraints_upper_bounds);
         user_model.number_jacobian_nonzeros = number_jacobian_nonzeros;
         user_model.jacobian_row_indices = std::move(jacobian_row_indices);
         user_model.jacobian_column_indices = std::move(jacobian_column_indices);
         user_model.jacobian = std::move(jacobian);
         user_model.initial_dual_iterate.resize(static_cast<size_t>(number_constraints), 0.);
      })

      .def("set_constraints_lower_bounds", [](PythonUserModel& user_model, std::vector<double> constraints_lower_bounds) {
         user_model.constraints_lower_bounds = std::move(constraints_lower_bounds);
      })

      .def("set_constraints_upper_bounds", [](PythonUserModel& user_model, std::vector<double> constraints_upper_bounds) {
         user_model.constraints_upper_bounds = std::move(constraints_upper_bounds);
      })

      .def("set_constraint_lower_bound", [](PythonUserModel& user_model, uno_int constraint_index, double lower_bound) {
         if (constraint_index < 0 || constraint_index >= user_model.number_constraints) {
            throw std::invalid_argument("Please specify a valid index.");
         }
         user_model.constraints_lower_bounds[static_cast<size_t>(constraint_index)] = lower_bound;
      })

      .def("set_constraint_upper_bound", [](PythonUserModel& user_model, uno_int constraint_index, double upper_bound) {
         if (constraint_index < 0 || constraint_index >= user_model.number_variables) {
            throw std::invalid_argument("Please specify a valid index.");
         }
         user_model.constraints_upper_bounds[static_cast<size_t>(constraint_index)] = upper_bound;
      })

      .def("set_lagrangian_hessian", [](PythonUserModel& user_model, uno_int number_hessian_nonzeros, char hessian_triangular_part,
            std::vector<uno_int> hessian_row_indices, std::vector<uno_int> hessian_column_indices, Hessian lagrangian_hessian) {
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
            if (hessian_row_indices.size() != static_cast<size_t>(number_hessian_nonzeros)) {
               throw std::invalid_argument("Dimension mismatch in hessian_row_indices.");
            }
            if (hessian_column_indices.size() != static_cast<size_t>(number_hessian_nonzeros)) {
               throw std::invalid_argument("Dimension mismatch in hessian_column_indices.");
            }
            // we only maintain the lower triangle
            if (hessian_triangular_part == UNO_LOWER_TRIANGLE) {
               user_model.hessian_row_indices = std::move(hessian_row_indices);
               user_model.hessian_column_indices = std::move(hessian_column_indices);
            }
            else {
               user_model.hessian_row_indices = std::move(hessian_column_indices);
               user_model.hessian_column_indices = std::move(hessian_row_indices);
            }
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

      .def("set_initial_primal_iterate", [](PythonUserModel& user_model, std::vector<double> initial_primal_iterate) {
         if (initial_primal_iterate.size() != static_cast<size_t>(user_model.number_variables)) {
            throw std::invalid_argument("Dimension mismatch in initial_primal_iterate.");
         }
         user_model.initial_primal_iterate = std::move(initial_primal_iterate);
      })

      .def("set_initial_dual_iterate", [](PythonUserModel& user_model, std::vector<double> initial_dual_iterate) {
         if (initial_dual_iterate.size() != static_cast<size_t>(user_model.number_constraints)) {
            throw std::invalid_argument("Dimension mismatch in initial_dual_iterate.");
         }
         user_model.initial_dual_iterate = std::move(initial_dual_iterate);
      });
   }
} // namespace