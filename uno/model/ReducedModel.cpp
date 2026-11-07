// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <exception>
#include <limits>
#include <stdexcept>
#include "ReducedModel.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"

namespace uno {
   namespace {
      constexpr size_t REMOVED = std::numeric_limits<size_t>::max();

      // original index -> reduced index (REMOVED if the element was removed)
      std::vector<size_t> inverse_map(const std::vector<size_t>& kept_indices, size_t original_size) {
         std::vector<size_t> reduced_indices(original_size, REMOVED);
         for (size_t reduced_index: Range(kept_indices.size())) {
            reduced_indices[kept_indices[reduced_index]] = reduced_index;
         }
         return reduced_indices;
      }

      // destination[k] = source[indices[k]]
      template <typename Source, typename Destination>
      void gather(const Source& source, const std::vector<size_t>& indices, Destination& destination) {
         for (size_t index: Range(indices.size())) {
            destination[index] = source[indices[index]];
         }
      }

      // destination[indices[k]] = source[k]
      template <typename Source, typename Destination>
      void scatter(const Source& source, const std::vector<size_t>& indices, Destination& destination) {
         for (size_t index: Range(indices.size())) {
            destination[indices[index]] = source[index];
         }
      }
   } // namespace

   // the original primals, with the fixed variables at their values
   template <typename Array>
   const Vector<double>& ReducedModel::expand_primals(const Array& x) const {
      scatter(x, this->reduction.free_variables, this->original_primals);
      return this->original_primals;
   }

   ReducedModel::ReducedModel(const Model& original_model, ModelReduction&& reduction):
         Model(original_model.name + " -> reduced", reduction.free_variables.size(), reduction.kept_constraints.size(),
            original_model.optimization_sense, original_model.lagrangian_sign_convention, original_model.base_indexing),
         model(original_model),
         reduction(std::move(reduction)),
         variables_lower_bounds(this->number_variables), variables_upper_bounds(this->number_variables),
         constraints_lower_bounds(this->number_constraints), constraints_upper_bounds(this->number_constraints),
         original_primals(original_model.number_variables), original_multipliers(original_model.number_constraints),
         original_constraints(original_model.number_constraints), original_gradient(original_model.number_variables),
         original_vector(original_model.number_variables), original_result(original_model.number_variables),
         original_jacobian_values(original_model.number_jacobian_nonzeros()) {
      const auto& free_variables = this->reduction.free_variables;
      const auto& kept_constraints = this->reduction.kept_constraints;
      // bounds
      gather(original_model.get_variables_lower_bounds(), free_variables, this->variables_lower_bounds);
      gather(original_model.get_variables_upper_bounds(), free_variables, this->variables_upper_bounds);
      gather(original_model.get_constraints_lower_bounds(), kept_constraints, this->constraints_lower_bounds);
      gather(original_model.get_constraints_upper_bounds(), kept_constraints, this->constraints_upper_bounds);
      for (size_t variable_index: this->reduction.fixed_variables) {
         this->original_primals[variable_index] = original_model.get_variables_lower_bounds()[variable_index];
      }
      this->partition_constraints(this->equality_constraints, this->inequality_constraints);

      // linear and nonlinear constraints
      const std::vector<size_t> reduced_constraint = inverse_map(kept_constraints, original_model.number_constraints);
      for (size_t constraint_index: original_model.get_linear_constraints()) {
         if (reduced_constraint[constraint_index] != REMOVED) {
            this->linear_constraints.push_back(reduced_constraint[constraint_index]);
         }
      }
      for (size_t constraint_index: original_model.get_nonlinear_constraints()) {
         if (reduced_constraint[constraint_index] != REMOVED) {
            this->nonlinear_constraints.push_back(reduced_constraint[constraint_index]);
         }
      }

      // Jacobian sparsity: the nonzeros of the kept constraints w.r.t. the free variables
      const std::vector<size_t> reduced_variable = inverse_map(free_variables, original_model.number_variables);
      const auto row_indices = original_model.get_jacobian_row_indices();
      const auto column_indices = original_model.get_jacobian_column_indices();
      for (size_t nonzero_index: Range(original_model.number_jacobian_nonzeros())) {
         const size_t row = reduced_constraint[static_cast<size_t>(row_indices[nonzero_index])];
         const size_t column = reduced_variable[static_cast<size_t>(column_indices[nonzero_index])];
         if (row != REMOVED && column != REMOVED) {
            this->jacobian_nonzeros.push_back(nonzero_index);
            this->jacobian_row_indices.push_back(static_cast<uno_int>(row));
            this->jacobian_column_indices.push_back(static_cast<uno_int>(column));
         }
      }

      // Hessian sparsity: the nonzeros w.r.t. the free variables
      if (original_model.has_hessian_matrix()) {
         const size_t number_hessian_nonzeros = original_model.number_hessian_nonzeros();
         std::vector<uno_int> hessian_rows(number_hessian_nonzeros), hessian_columns(number_hessian_nonzeros);
         original_model.compute_hessian_sparsity(view(hessian_rows.data(), number_hessian_nonzeros),
            view(hessian_columns.data(), number_hessian_nonzeros), 0);
         for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
            const size_t row = reduced_variable[static_cast<size_t>(hessian_rows[nonzero_index])];
            const size_t column = reduced_variable[static_cast<size_t>(hessian_columns[nonzero_index])];
            if (row != REMOVED && column != REMOVED) {
               this->hessian_nonzeros.push_back(nonzero_index);
               this->hessian_row_indices.push_back(static_cast<uno_int>(row));
               this->hessian_column_indices.push_back(static_cast<uno_int>(column));
            }
         }
         this->original_hessian_values.resize(number_hessian_nonzeros);
      }
      DISCRETE << "Preprocessing: " << this->reduction.fixed_variables.size() << " fixed variable(s) eliminated, " <<
         this->reduction.removed_constraints.size() << " constraint(s) of fixed variables only removed\n";
   }

   ProblemType ReducedModel::get_problem_type() const {
      return this->model.get_problem_type();
   }

   bool ReducedModel::has_jacobian_operator() const {
      return false;
   }

   bool ReducedModel::has_jacobian_transposed_operator() const {
      return false;
   }

   bool ReducedModel::has_hessian_operator() const {
      return this->model.has_hessian_operator();
   }

   bool ReducedModel::has_hessian_matrix() const {
      return this->model.has_hessian_matrix();
   }

   double ReducedModel::evaluate_objective(const Vector<double>& x) const {
      return this->model.evaluate_objective(this->expand_primals(x));
   }

   void ReducedModel::evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const {
      this->model.evaluate_constraints(this->expand_primals(x), this->original_constraints);
      gather(this->original_constraints, this->reduction.kept_constraints, constraints);
   }

   void ReducedModel::evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const {
      this->original_gradient.fill(0.);
      this->model.evaluate_objective_gradient(this->expand_primals(x), this->original_gradient);
      gather(this->original_gradient, this->reduction.free_variables, gradient);
   }

   View<const uno_int> ReducedModel::get_jacobian_row_indices() const {
      return view(this->jacobian_row_indices.data(), this->jacobian_row_indices.size());
   }

   View<const uno_int> ReducedModel::get_jacobian_column_indices() const {
      return view(this->jacobian_column_indices.data(), this->jacobian_column_indices.size());
   }

   void ReducedModel::compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const {
      for (size_t nonzero_index: Range(this->hessian_nonzeros.size())) {
         row_indices[nonzero_index] = this->hessian_row_indices[nonzero_index] + solver_indexing;
         column_indices[nonzero_index] = this->hessian_column_indices[nonzero_index] + solver_indexing;
      }
   }

   void ReducedModel::evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const {
      this->model.evaluate_jacobian(this->expand_primals(x), this->original_jacobian_values.data());
      gather(this->original_jacobian_values, this->jacobian_nonzeros, jacobian_values);
   }

   // the multipliers of the bound constraints are 0 (their Hessian vanishes anyway)
   void ReducedModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         View<double> hessian_values) const {
      scatter(multipliers, this->reduction.kept_constraints, this->original_multipliers);
      this->model.evaluate_lagrangian_hessian(this->expand_primals(x), objective_multiplier, this->original_multipliers,
         this->original_hessian_values.view());
      gather(this->original_hessian_values, this->hessian_nonzeros, hessian_values);
   }

   void ReducedModel::compute_jacobian_vector_product(const double* /*x*/, const double* /*vector*/, double* /*result*/) const {
      throw std::runtime_error("ReducedModel::compute_jacobian_vector_product not implemented");
   }

   void ReducedModel::compute_jacobian_transposed_vector_product(const double* /*x*/, const double* /*vector*/, double* /*result*/) const {
      throw std::runtime_error("ReducedModel::compute_jacobian_transposed_vector_product not implemented");
   }

   // the fixed components of the vector are 0
   void ReducedModel::compute_hessian_vector_product(View<const double> x, View<const double> vector, double objective_multiplier,
         const Vector<double>& multipliers, View<double> result) const {
      scatter(vector, this->reduction.free_variables, this->original_vector);
      scatter(multipliers, this->reduction.kept_constraints, this->original_multipliers);
      this->model.compute_hessian_vector_product(this->expand_primals(x).view(), this->original_vector.view(),
         objective_multiplier, this->original_multipliers, this->original_result.view());
      gather(this->original_result, this->reduction.free_variables, result);
   }

   const std::vector<double>& ReducedModel::get_variables_lower_bounds() const {
      return this->variables_lower_bounds;
   }

   const std::vector<double>& ReducedModel::get_variables_upper_bounds() const {
      return this->variables_upper_bounds;
   }

   const std::vector<double>& ReducedModel::get_constraints_lower_bounds() const {
      return this->constraints_lower_bounds;
   }

   const std::vector<double>& ReducedModel::get_constraints_upper_bounds() const {
      return this->constraints_upper_bounds;
   }

   const std::vector<size_t>& ReducedModel::get_equality_constraints() const {
      return this->equality_constraints;
   }

   const std::vector<size_t>& ReducedModel::get_inequality_constraints() const {
      return this->inequality_constraints;
   }

   const Collection<size_t>& ReducedModel::get_linear_constraints() const {
      return this->linear_constraints_collection;
   }

   const Collection<size_t>& ReducedModel::get_nonlinear_constraints() const {
      return this->nonlinear_constraints_collection;
   }

   void ReducedModel::initial_primal_point(Vector<double>& x) const {
      Vector<double> original_x(this->model.number_variables);
      this->model.initial_primal_point(original_x);
      gather(original_x, this->reduction.free_variables, x);
   }

   void ReducedModel::initial_dual_point(Vector<double>& multipliers) const {
      Vector<double> original_initial_multipliers(this->model.number_constraints);
      this->model.initial_dual_point(original_initial_multipliers);
      gather(original_initial_multipliers, this->reduction.kept_constraints, multipliers);
   }

   void ReducedModel::postprocess_solution(Iterate& iterate, Evaluations& evaluations) const {
      const Iterate reduced_iterate = iterate;
      iterate.set_dimensions(this->model.number_variables, this->model.number_constraints);
      iterate.primals = this->expand_primals(reduced_iterate.primals);

      // the cached evaluations live in the reduced space; the objective value is unchanged
      const double objective = evaluations.objective;
      const bool is_objective_computed = evaluations.is_objective_computed;
      evaluations = Evaluations(this->model);
      evaluations.objective = objective;
      evaluations.is_objective_computed = is_objective_computed;

      this->recover_multipliers(reduced_iterate, iterate, evaluations);
      this->model.postprocess_solution(iterate, evaluations);
   }

   size_t ReducedModel::number_jacobian_nonzeros() const {
      return this->jacobian_nonzeros.size();
   }

   size_t ReducedModel::number_hessian_nonzeros() const {
      return this->hessian_nonzeros.size();
   }

   size_t ReducedModel::number_model_objective_evaluations() const {
      return this->model.number_model_objective_evaluations();
   }

   size_t ReducedModel::number_model_constraints_evaluations() const {
      return this->model.number_model_constraints_evaluations();
   }

   size_t ReducedModel::number_model_objective_gradient_evaluations() const {
      return this->model.number_model_objective_gradient_evaluations();
   }

   size_t ReducedModel::number_model_jacobian_evaluations() const {
      return this->model.number_model_jacobian_evaluations();
   }

   size_t ReducedModel::number_model_hessian_evaluations() const {
      return this->model.number_model_hessian_evaluations();
   }

   void ReducedModel::reset_number_evaluations() const {
      this->model.reset_number_evaluations();
   }

   // with the Lagrangian ρ f(x) - y^T c(x) - z^T x, a fixed variable x_i gets the bound multiplier
   // z_i = [ρ ∇f(x) - ∇c(x) y]_i, split into its lower (z_i > 0) and upper (z_i < 0) parts. The removed constraints get 0
   void ReducedModel::recover_multipliers(const Iterate& reduced_iterate, Iterate& iterate, Evaluations& evaluations) const {
      Multipliers& multipliers = iterate.multipliers;
      multipliers.reset();
      scatter(reduced_iterate.multipliers.constraints, this->reduction.kept_constraints, multipliers.constraints);
      scatter(reduced_iterate.multipliers.lower_bounds, this->reduction.free_variables, multipliers.lower_bounds);
      scatter(reduced_iterate.multipliers.upper_bounds, this->reduction.free_variables, multipliers.upper_bounds);
      try {
         // the bound multipliers of the fixed variables are still 0
         Vector<double> lagrangian_gradient(this->model.number_variables);
         this->model.evaluate_lagrangian_gradient(iterate.primals, multipliers.constraints, multipliers.lower_bounds,
            multipliers.upper_bounds, iterate.objective_multiplier, evaluations, lagrangian_gradient.view());
         for (size_t variable_index: this->reduction.fixed_variables) {
            multipliers.lower_bounds[variable_index] = std::max(lagrangian_gradient[variable_index], 0.);
            multipliers.upper_bounds[variable_index] = std::min(lagrangian_gradient[variable_index], 0.);
         }
      }
      catch (const std::exception& exception) {
         WARNING << "Preprocessing: the multipliers of the fixed variables could not be computed (" << exception.what() << ")\n";
      }
   }
} // namespace
