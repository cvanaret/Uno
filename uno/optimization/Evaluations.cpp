// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Evaluations.hpp"
#include "model/Model.hpp"

namespace uno {
   Evaluations::Evaluations(size_t number_variables, size_t number_constraints):
      constraints(number_constraints),
      objective_gradient(number_variables) {
   }

   void Evaluations::evaluate_objective(const Model& model, const Vector<double>& primals) {
      this->objective = model.evaluate_objective(primals);
   }

   void Evaluations::evaluate_constraints(const Model& model, const Vector<double>& primals) {
      model.evaluate_constraints(primals, this->constraints);
   }

   void Evaluations::evaluate_objective_gradient(const Model& model, const Vector<double>& primals) {
      model.evaluate_objective_gradient(primals, this->objective_gradient);
   }

   void Evaluations::compute_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const {

   }

   void Evaluations::compute_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const {

   }
} // namespace