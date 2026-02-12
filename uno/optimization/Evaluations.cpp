// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include "Evaluations.hpp"
#include "linear_algebra/COOSparsity.hpp"
#include "model/Model.hpp"
#include "optimization/EvaluationErrors.hpp"

namespace uno {
   Evaluations::Evaluations(const Model& model, const COOSparsity* jacobian_sparsity):
      constraints(model.number_constraints),
      objective_gradient(model.number_variables),
      jacobian_values(model.number_jacobian_nonzeros()),
      jacobian_sparsity(jacobian_sparsity) {
   }

   void Evaluations::evaluate_objective(const Model& model, const Vector<double>& primals) {
      if (!this->is_objective_computed) {
         this->objective = model.evaluate_objective(primals);
         // check finiteness
         if (!is_finite(this->objective)) {
            throw FunctionEvaluationError();
         }
         this->is_objective_computed = true;
      }
   }

   void Evaluations::evaluate_constraints(const Model& model, const Vector<double>& primals) {
      if (!this->are_constraints_computed) {
         if (model.is_constrained()) {
            model.evaluate_constraints(primals, this->constraints);
            // check finiteness
            if (std::any_of(this->constraints.begin(), this->constraints.end(), [](double constraint_value) {
               return !is_finite(constraint_value);
            })) {
               throw FunctionEvaluationError();
            }
         }
         this->are_constraints_computed = true;
      }
   }

   void Evaluations::evaluate_objective_gradient(const Model& model, const Vector<double>& primals) {
      if (!this->is_objective_gradient_computed) {
         this->objective_gradient.fill(0.);
         model.evaluate_objective_gradient(primals, this->objective_gradient);
         // check finiteness
         if (std::any_of(this->objective_gradient.begin(), this->objective_gradient.end(), [](double derivative) {
            return !is_finite(derivative);
         })) {
            throw GradientEvaluationError();
         }
         this->is_objective_gradient_computed = true;
      }
   }

   void Evaluations::evaluate_jacobian(const Model& model, const Vector<double>& primals) {
      model.evaluate_jacobian(primals, this->jacobian_values.data());
      this->is_jacobian_computed = true;
   }

   void Evaluations::compute_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      result.fill(0.);
      const size_t number_jacobian_nonzeros = this->jacobian_sparsity->row_indices.size();
      for (size_t nonzero_index: Range(number_jacobian_nonzeros)) {
         const size_t constraint_index = static_cast<size_t>(this->jacobian_sparsity->row_indices[nonzero_index]);
         const size_t variable_index = static_cast<size_t>(this->jacobian_sparsity->column_indices[nonzero_index]);
         const double derivative = this->jacobian_values[nonzero_index];
         assert(variable_index < vector.size());
         assert(constraint_index < result.size());

         result[constraint_index] += derivative * vector[variable_index];
      }
   }

   void Evaluations::compute_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      result.fill(0.);
      const size_t number_jacobian_nonzeros = this->jacobian_sparsity->row_indices.size();
      for (size_t nonzero_index: Range(number_jacobian_nonzeros)) {
         const size_t constraint_index = static_cast<size_t>(this->jacobian_sparsity->row_indices[nonzero_index]);
         const size_t variable_index = static_cast<size_t>(this->jacobian_sparsity->column_indices[nonzero_index]);
         const double derivative = this->jacobian_values[nonzero_index];
         assert(constraint_index < vector.size());
         assert(variable_index < result.size());

         result[variable_index] += derivative * vector[constraint_index];
      }
   }
} // namespace