// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ScaledModel.hpp"
#include "Model.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   ScaledModel::ScaledModel(const Model& original_model, const Vector<double>& initial_primals, const Options& options):
         Model(original_model.name + " -> scaled", original_model.number_variables, original_model.number_constraints,
               original_model.optimization_sense, original_model.lagrangian_sign_convention, original_model.base_indexing),
         model(original_model),
         scaling(this->model.number_constraints, options.get_double("function_scaling_threshold")),
         constraints_lower_bounds(this->model.get_constraints_lower_bounds()),
         constraints_upper_bounds(this->model.get_constraints_upper_bounds()) {
      if (options.get_bool("use_function_scaling")) {
         // evaluate the gradients at the current point
         Vector<double> objective_gradient(this->model.number_variables);
         Vector<double> jacobian_values(this->model.number_jacobian_nonzeros());
         try {
            this->model.evaluate_objective_gradient(initial_primals, objective_gradient);
            this->model.evaluate_jacobian(initial_primals, jacobian_values.data());
            this->scaling.compute(this->model, objective_gradient, jacobian_values);
         }
         catch (const EvaluationError&) {
            DEBUG << "The gradients could not be evaluated at the initial point, functions will not be scaled\n";
         }
      }

      // compute the lower and upper bounds of the constraints
      if (this->scaling.are_constraints_scaled()) {
         const auto& constraints_scaling = this->scaling.get_constraint_scaling();
         for (size_t constraint_index: Range(this->number_constraints)) {
            this->constraints_lower_bounds[constraint_index] *= constraints_scaling[constraint_index];
            this->constraints_upper_bounds[constraint_index] *= constraints_scaling[constraint_index];
         }

         // preallocate the temporary multipliers
         this->scaled_multipliers.resize(this->number_constraints);
      }
   }

   ProblemType ScaledModel::get_problem_type() const {
      return this->model.get_problem_type();
   }

   bool ScaledModel::has_jacobian_operator() const {
      return this->model.has_jacobian_operator();
   }

   bool ScaledModel::has_jacobian_transposed_operator() const {
      return this->model.has_jacobian_transposed_operator();
   }

   bool ScaledModel::has_hessian_operator() const {
      return this->model.has_hessian_operator();
   }

   bool ScaledModel::has_hessian_matrix() const {
      return this->model.has_hessian_matrix();
   }

   double ScaledModel::evaluate_objective(const Vector<double>& x) const {
      const double objective = this->model.evaluate_objective(x);
      return this->scaling.get_objective_scaling() * objective;
   }

   void ScaledModel::evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const {
      this->model.evaluate_objective_gradient(x, gradient);
      if (this->scaling.is_objective_scaled()) {
         gradient.scale(this->scaling.get_objective_scaling());
      }
   }

   View<const uno_int> ScaledModel::get_jacobian_row_indices() const {
      return this->model.get_jacobian_row_indices();
   }

   View<const uno_int> ScaledModel::get_jacobian_column_indices() const {
      return this->model.get_jacobian_column_indices();
   }

   void ScaledModel::compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const {
      this->model.compute_hessian_sparsity(row_indices, column_indices, solver_indexing);
   }

   void ScaledModel::evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const {
      this->model.evaluate_constraints(x, constraints);
      if (this->scaling.are_constraints_scaled()) {
         const auto& constraints_scaling = this->scaling.get_constraint_scaling();
         for (size_t constraint_index: Range(this->number_constraints)) {
            constraints[constraint_index] *= constraints_scaling[constraint_index];
         }
      }
   }

   void ScaledModel::evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const {
      this->model.evaluate_jacobian(x, jacobian_values);
      // scale each term of the Jacobian, depending on which row/constraint it belongs to
      if (this->scaling.are_constraints_scaled()) {
         const auto& jacobian_row_indices = this->model.get_jacobian_row_indices();
         const auto& constraints_scaling = this->scaling.get_constraint_scaling();
         for (size_t nonzero_index: Range(this->model.number_jacobian_nonzeros())) {
            const size_t constraint_index = static_cast<size_t>(jacobian_row_indices[nonzero_index]);
            jacobian_values[nonzero_index] *= constraints_scaling[constraint_index];
         }
      }
   }

   void ScaledModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         View<double> hessian_values) const {
      // scale the objective and constraint multipliers (in a lazy way for the multipliers)
      const double scaled_objective_multiplier = objective_multiplier * this->scaling.get_objective_scaling();
      if (this->scaling.are_constraints_scaled()) {
         const auto& constraints_scaling = this->scaling.get_constraint_scaling();
         for (size_t constraint_index: Range(this->number_constraints)) {
            this->scaled_multipliers[constraint_index] = constraints_scaling[constraint_index] * multipliers[constraint_index];
         }
         this->model.evaluate_lagrangian_hessian(x, scaled_objective_multiplier, this->scaled_multipliers, hessian_values);
      }
      else {
         this->model.evaluate_lagrangian_hessian(x, scaled_objective_multiplier, multipliers, hessian_values);
      }
   }

   void ScaledModel::compute_jacobian_vector_product(const double* /*x*/, const double* /*vector*/, double* /*result*/) const {
      throw std::runtime_error("ScaledModel::compute_jacobian_vector_product not implemented yet");
   }

   void ScaledModel::compute_jacobian_transposed_vector_product(const double* /*x*/, const double* /*vector*/, double* /*result*/) const {
      throw std::runtime_error("ScaledModel::compute_jacobian_transposed_vector_product not implemented yet");
   }

   void ScaledModel::compute_hessian_vector_product(View<const double> x, View<const double> vector, double objective_multiplier,
         const Vector<double>& multipliers, View<double> result) const {
      // scale the objective and constraint multipliers (in a lazy way for the multipliers)
      const double scaled_objective_multiplier = objective_multiplier * this->scaling.get_objective_scaling();
      if (this->scaling.are_constraints_scaled()) {
         const auto& constraints_scaling = this->scaling.get_constraint_scaling();
         for (size_t constraint_index: Range(this->number_constraints)) {
            this->scaled_multipliers[constraint_index] = constraints_scaling[constraint_index] * multipliers[constraint_index];
         }
         this->model.compute_hessian_vector_product(x, vector, scaled_objective_multiplier, this->scaled_multipliers, result);
      }
      else {
         this->model.compute_hessian_vector_product(x, vector, scaled_objective_multiplier, multipliers, result);
      }
   }

   const std::vector<double>& ScaledModel::get_variables_lower_bounds() const {
      return this->model.get_variables_lower_bounds();
   }

   const std::vector<double>& ScaledModel::get_variables_upper_bounds() const {
      return this->model.get_variables_upper_bounds();
   }

   const Vector<size_t>& ScaledModel::get_fixed_variables() const {
      return this->model.get_fixed_variables();
   }

   const std::vector<double>& ScaledModel::get_constraints_lower_bounds() const {
      return this->constraints_lower_bounds;
   }

   const std::vector<double>& ScaledModel::get_constraints_upper_bounds() const {
      return this->constraints_upper_bounds;
   }

   const Collection<size_t>& ScaledModel::get_equality_constraints() const {
      return this->model.get_equality_constraints();
   }

   const Collection<size_t>& ScaledModel::get_inequality_constraints() const {
      return this->model.get_inequality_constraints();
   }

   const Collection<size_t>& ScaledModel::get_linear_constraints() const {
      return this->model.get_linear_constraints();
   }

   const Collection<size_t>& ScaledModel::get_nonlinear_constraints() const {
      return this->model.get_nonlinear_constraints();
   }

   void ScaledModel::initial_primal_point(Vector<double>& x) const {
      this->model.initial_primal_point(x);
   }
   
   void ScaledModel::initial_dual_point(Vector<double>& multipliers) const {
      this->model.initial_dual_point(multipliers);
   }

   void ScaledModel::postprocess_solution(Iterate& iterate, Evaluations& evaluations) const {
      // unscale the constraints and the constraint multipliers
      if (scaling.is_objective_scaled() || scaling.are_constraints_scaled()) {
         const double objective_scaling = this->scaling.get_objective_scaling();
         const auto& constraints_scaling = this->scaling.get_constraint_scaling();
         for (size_t constraint_index: Range(this->model.number_constraints)) {
            evaluations.constraints[constraint_index] /= constraints_scaling[constraint_index];
            iterate.multipliers.constraints[constraint_index] *= constraints_scaling[constraint_index] / objective_scaling;
         }
      }


      // unscale the objective value and the bound multipliers
      if (scaling.is_objective_scaled()) {
         const double objective_scaling = this->scaling.get_objective_scaling();
         if (evaluations.is_objective_computed) {
            evaluations.objective /= objective_scaling;
         }

         for (size_t variable_index: Range(this->model.number_variables)) {
            iterate.multipliers.lower_bounds[variable_index] /= objective_scaling;
            iterate.multipliers.upper_bounds[variable_index] /= objective_scaling;
         }
      }

      this->model.postprocess_solution(iterate, evaluations);
   }

   size_t ScaledModel::number_jacobian_nonzeros() const {
      return this->model.number_jacobian_nonzeros();
   }

   size_t ScaledModel::number_hessian_nonzeros() const {
      return this->model.number_hessian_nonzeros();
   }

   size_t ScaledModel::number_model_objective_evaluations() const {
      return this->model.number_model_objective_evaluations();
   }

   size_t ScaledModel::number_model_constraints_evaluations() const {
      return this->model.number_model_constraints_evaluations();
   }

   size_t ScaledModel::number_model_objective_gradient_evaluations() const {
      return this->model.number_model_objective_gradient_evaluations();
   }

   size_t ScaledModel::number_model_jacobian_evaluations() const {
      return this->model.number_model_jacobian_evaluations();
   }

   size_t ScaledModel::number_model_hessian_evaluations() const {
      return this->model.number_model_hessian_evaluations();
   }

   void ScaledModel::reset_number_evaluations() const {
      this->model.reset_number_evaluations();
   }
} // namespace