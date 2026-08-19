// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ExactHessian.hpp"
#include "model/Model.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "tools/Logger.hpp"

namespace uno {
   // Hessian buffer
   HessianBuffer::HessianBuffer(const Model& model): hessian_buffer(model.number_hessian_nonzeros()) {
   }

   void HessianBuffer::evaluate_hessian(const Model& model, const Vector<double>& primal_variables, double objective_multiplier,
         const Vector<double>& constraint_multipliers, View<double> hessian_values) {
      assert(this->hessian_buffer.size() == hessian_values.size());

      // try to evaluate the Hessian. Upon failure, keep the previous one
      try {
         model.evaluate_lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, this->hessian_buffer.view());
         // success: copy the new Hessian into the Hessian values
         hessian_values = this->hessian_buffer;
      }
      catch (const HessianEvaluationError&) {
         if (this->first_evaluation) {
            hessian_values.fill(1.);
            DEBUG << "The initial Hessian could not be evaluated, using the identity\n";
         }
         else {
            // keep the existing Hessian
            DEBUG << "The Hessian could not be evaluated, using the previous one\n";
         }
      }
      this->first_evaluation = false;
   }

   // exact Hessian

   ExactHessian::ExactHessian(const Model& model): HessianModel("exact"), model(model), hessian_buffer(model) {
   }

   bool ExactHessian::has_hessian_operator() const {
      return this->model.has_hessian_operator();
   }

   bool ExactHessian::has_hessian_matrix() const {
      return this->model.has_hessian_matrix();
   }

   bool ExactHessian::has_curvature() const {
      return (this->model.get_problem_type() != ProblemType::LINEAR);
   }

   size_t ExactHessian::number_nonzeros() const {
      return this->model.number_hessian_nonzeros();
   }

   void ExactHessian::compute_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing) const {
      // Hessian sparsity of the model
      this->model.compute_hessian_sparsity(row_indices, column_indices, solver_indexing);
   }

   bool ExactHessian::is_positive_definite() const {
      return false;
   }

   void ExactHessian::initialize_statistics(Statistics& /*statistics*/) const {
   }

   void ExactHessian::notify_trial_iterate(Statistics& /*statistics*/, const Iterate& /*current_iterate*/,
         const Iterate& /*trial_iterate*/, Evaluations& /*current_evaluations*/, Evaluations& /*trial_evaluations*/) {
   }

   void ExactHessian::evaluate_hessian(Statistics& /*statistics*/, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, View<double> hessian_values) {
      this->model.evaluate_lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, hessian_values);
      ++this->evaluation_count;
   }

   void ExactHessian::compute_hessian_vector_product(const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* result) {
      this->model.compute_hessian_vector_product(x, vector, objective_multiplier, constraint_multipliers, result);
      ++this->evaluation_count;
   }
} // namespace