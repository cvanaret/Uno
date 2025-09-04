// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include "../unopy.hpp"
#include "symbolic/Range.hpp"
#include "tools/PointerWrapper.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(uno::PointerWrapper<const uno::Vector<double>>);

namespace uno {
   using PyObjective = std::function<double(PointerWrapper<const Vector<double>>)>;
   using PyObjectiveGradient = std::function<void(PointerWrapper<const Vector<double>>, PointerWrapper<Vector<double>>)>;

   void define_Model(py::module& module) {
      py::class_<PythonUserModel>(module, "Model")
      // constructor
      // https://stackoverflow.com/a/62310838/16037994
      .def(py::init<>([](char problem_type, int32_t number_variables, std::vector<double>& variables_lower_bounds,
         std::vector<double>& variables_upper_bounds, int32_t base_indexing) {
         return PythonUserModel(problem_type, number_variables, variables_lower_bounds.data(), variables_upper_bounds.data(),
            base_indexing);
      }), "Constructor")

      // methods
      .def("set_objective", [](PythonUserModel& user_model, int32_t optimization_sense, const PyObjective& objective_function,
            const PyObjectiveGradient& objective_gradient) {
         if (optimization_sense != UNO_MINIMIZE && optimization_sense != UNO_MAXIMIZE) {
            throw std::runtime_error("Please specify a valid objective sense.");
         }
         user_model.optimization_sense = optimization_sense;
         //user_model.objective_function = objective_function;
         //user_model.objective_gradient = objective_gradient;
      })

      .def("set_constraints", [](PythonUserModel& user_model, int32_t number_constraints, Constraints&& constraint_functions,
            std::vector<double>& constraints_lower_bounds, std::vector<double>& constraints_upper_bounds, int32_t number_jacobian_nonzeros,
            const std::vector<int32_t>& jacobian_row_indices, const std::vector<int32_t>& jacobian_column_indices, Jacobian&& constraint_jacobian) {
         if (number_constraints <= 0) {
            throw std::runtime_error("Please specify a positive number of constraints.");
         }
         user_model.number_constraints = number_constraints;
         user_model.constraint_functions = constraint_functions;
         user_model.constraints_lower_bounds = constraints_lower_bounds.data();
         user_model.constraints_upper_bounds = constraints_upper_bounds.data();
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

      .def("set_initial_primal_iterate", [](PythonUserModel& user_model, std::vector<double>& initial_primal_iterate) {
         user_model.initial_primal_iterate = initial_primal_iterate;
      });
   }
} // namespace