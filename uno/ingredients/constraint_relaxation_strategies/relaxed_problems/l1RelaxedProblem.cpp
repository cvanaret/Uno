// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "l1RelaxedProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "optimization/EvaluationSpace.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   l1RelaxedProblem::l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient,
         bool relax_linear_constraints):
      // call delegating constructor
      l1RelaxedProblem(model, ElasticVariables::generate(model, relax_linear_constraints), objective_multiplier,
         constraint_violation_coefficient) {
   }

   // delegating constructor
   l1RelaxedProblem::l1RelaxedProblem(const Model& model, ElasticVariables&& elastic_variables, double objective_multiplier,
            double constraint_violation_coefficient):
         OptimizationProblem(model, model.number_variables + elastic_variables.size(), model.number_constraints),
         elastic_variables(elastic_variables),
         number_elastic_variables(elastic_variables.size()),
         objective_multiplier(objective_multiplier),
         constraint_violation_coefficient(constraint_violation_coefficient) {
   }

   double l1RelaxedProblem::get_objective_multiplier() const {
      return this->objective_multiplier;
   }

   void l1RelaxedProblem::set_proximal_coefficient(double proximal_coefficient) {
      this->proximal_coefficient = proximal_coefficient;
   }

   void l1RelaxedProblem::set_proximal_center(double* proximal_center) {
      this->proximal_center = proximal_center;
   }

   void l1RelaxedProblem::evaluate_constraints(Iterate& iterate, Vector<double>& constraints) const {
      iterate.evaluate_constraints(this->model);
      constraints = iterate.evaluations.constraints;

      // add the contribution of the elastic variables
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
         constraints[constraint_index] -= iterate.primals[elastic_index];
      }
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
         constraints[constraint_index] += iterate.primals[elastic_index];
      }
   }

   void l1RelaxedProblem::evaluate_objective_gradient(Iterate& iterate, double* objective_gradient) const {
      // scale nabla f(x) by rho
      if (this->objective_multiplier != 0.) {
         iterate.evaluate_objective_gradient(this->model);
         // TODO change this
         for (size_t index: Range(this->model.number_variables)) {
            objective_gradient[index] = this->objective_multiplier * iterate.evaluations.objective_gradient[index];
         }
      }
      else {
         for (size_t index: Range(this->model.number_variables)) {
            objective_gradient[index] = 0.;
         }
      }

      // constraint violation (through elastic variables) contribution
      for (size_t elastic_index: Range(this->model.number_variables, this->number_variables)) {
         objective_gradient[elastic_index] = this->constraint_violation_coefficient;
      }

      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double proximal_term = this->proximal_coefficient * scaling * scaling * (iterate.primals[variable_index] - this->proximal_center[variable_index]);
            objective_gradient[variable_index] += proximal_term;
         }
      }
   }

   size_t l1RelaxedProblem::number_jacobian_nonzeros() const {
      return this->model.number_jacobian_nonzeros() + this->number_elastic_variables;
   }

   bool l1RelaxedProblem::has_curvature(const HessianModel& hessian_model) const {
      // the l1 relaxation does not introduce curvature
      return hessian_model.has_curvature();
   }

   size_t l1RelaxedProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      size_t number_nonzeros = hessian_model.number_nonzeros();
      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         number_nonzeros += this->model.number_variables;
      }
      return number_nonzeros;
   }

   void l1RelaxedProblem::compute_constraint_jacobian_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing,
         MatrixOrder matrix_order) const {
      this->model.compute_constraint_jacobian_sparsity(row_indices, column_indices, solver_indexing, matrix_order);

      // add the contribution of the elastic variables
      size_t nonzero_index = this->model.number_jacobian_nonzeros();
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
         row_indices[nonzero_index] = static_cast<int>(constraint_index) + solver_indexing;
         column_indices[nonzero_index] = static_cast<uno_int>(elastic_index) + solver_indexing;
         ++nonzero_index;
      }
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
         row_indices[nonzero_index] = static_cast<int>(constraint_index) + solver_indexing;
         column_indices[nonzero_index] = static_cast<uno_int>(elastic_index) + solver_indexing;
         ++nonzero_index;
      }
   }

   void l1RelaxedProblem::compute_hessian_sparsity(const HessianModel& hessian_model, uno_int* row_indices,
         uno_int* column_indices, uno_int solver_indexing) const {
      hessian_model.compute_sparsity(row_indices, column_indices, solver_indexing);

      // diagonal proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         size_t current_index = hessian_model.number_nonzeros();
         for (size_t variable_index: Range(this->model.number_variables)) {
            row_indices[current_index] = static_cast<int>(variable_index) + solver_indexing;
            column_indices[current_index] = static_cast<int>(variable_index) + solver_indexing;
            ++current_index;
         }
      }
   }

   void l1RelaxedProblem::evaluate_constraint_jacobian(Iterate& iterate, double* jacobian_values) const {
      this->model.evaluate_constraint_jacobian(iterate.primals, jacobian_values);

      // add the contribution of the elastic variables
      size_t nonzero_index = this->model.number_jacobian_nonzeros();
      for ([[maybe_unused]] const auto _: this->elastic_variables.positive) {
         jacobian_values[nonzero_index] = -1.;
         ++nonzero_index;
      }
      for ([[maybe_unused]] const auto _: this->elastic_variables.negative) {
         jacobian_values[nonzero_index] = 1.;
         ++nonzero_index;
      }
   }

   // Lagrangian gradient split in two parts: objective contribution and constraints' contribution
   void l1RelaxedProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient,
         const EvaluationSpace& evaluation_space, Iterate& iterate) const {
      lagrangian_gradient.objective_contribution.fill(0.);
      lagrangian_gradient.constraints_contribution.fill(0.);

      // objective gradient
      lagrangian_gradient.objective_contribution = iterate.evaluations.objective_gradient;

      // ∇c(x_k) λ_k
      evaluation_space.compute_constraint_jacobian_transposed_vector_product(iterate.multipliers.constraints,
         lagrangian_gradient.constraints_contribution);
      lagrangian_gradient.constraints_contribution = -lagrangian_gradient.constraints_contribution;

      // bound constraints of original variables
      for (size_t variable_index: Range(this->model.number_variables)) {
         lagrangian_gradient.constraints_contribution[variable_index] -= (iterate.multipliers.lower_bounds[variable_index] +
            iterate.multipliers.upper_bounds[variable_index]);
      }

      // elastic variables
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
         lagrangian_gradient.constraints_contribution[elastic_index] += this->constraint_violation_coefficient +
            iterate.multipliers.constraints[constraint_index] - iterate.multipliers.lower_bounds[elastic_index];
      }
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
         lagrangian_gradient.constraints_contribution[elastic_index] += this->constraint_violation_coefficient -
            iterate.multipliers.constraints[constraint_index] - iterate.multipliers.lower_bounds[elastic_index];
      }

      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double proximal_term = this->proximal_coefficient * scaling * scaling;
            lagrangian_gradient.constraints_contribution[variable_index] += proximal_term * (iterate.primals[variable_index] -
               this->proximal_center[variable_index]);
         }
      }
   }

   void l1RelaxedProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, double* hessian_values) const {
      hessian_model.evaluate_hessian(statistics, primal_variables, this->get_objective_multiplier(),
         multipliers.constraints, hessian_values);

      // proximal contribution
      size_t nonzero_index = hessian_model.number_nonzeros();
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double proximal_term = this->proximal_coefficient * scaling * scaling;
            hessian_values[nonzero_index] = proximal_term;
            ++nonzero_index;
         }
      }
   }

   void l1RelaxedProblem::compute_hessian_vector_product(HessianModel& hessian_model, const double* x, const double* vector,
         const Multipliers& multipliers, double* result) const {
      hessian_model.compute_hessian_vector_product(x, vector, this->get_objective_multiplier(), multipliers.constraints, result);

      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double proximal_term = this->proximal_coefficient * scaling * scaling;
            result[variable_index] += proximal_term * vector[variable_index];
         }
      }
   }

   SolutionStatus l1RelaxedProblem::check_first_order_convergence(const Iterate& current_iterate, double primal_tolerance,
         double dual_tolerance) const {
      // evaluate termination conditions based on optimality conditions
      const bool feasibility_stationarity = (current_iterate.residuals.stationarity <= dual_tolerance);
      const bool primal_feasibility = (current_iterate.primal_feasibility <= primal_tolerance);
      const bool feasibility_complementarity = (current_iterate.residuals.complementarity <= dual_tolerance);
      const bool no_trivial_duals = current_iterate.multipliers.not_all_zero(this->model.number_variables, dual_tolerance);

      DEBUG << "\nTermination criteria for primal-dual tolerances = (" << primal_tolerance << ", " << dual_tolerance << "):\n";
      DEBUG << "Primal feasibility: " << std::boolalpha << primal_feasibility << '\n';
      DEBUG << "Feasibility stationarity: " << std::boolalpha << feasibility_stationarity << '\n';
      DEBUG << "Feasibility complementarity: " << std::boolalpha << feasibility_complementarity << '\n';
      DEBUG << "Not all zero multipliers: " << std::boolalpha << no_trivial_duals << "\n\n";

      if (this->model.is_constrained() && feasibility_stationarity && !primal_feasibility && feasibility_complementarity &&
            no_trivial_duals) {
         // no primal feasibility, stationary point of constraint violation
         return SolutionStatus::INFEASIBLE_STATIONARY_POINT;
      }
      return SolutionStatus::NOT_OPTIMAL;
   }

   double l1RelaxedProblem::variable_lower_bound(size_t variable_index) const {
      if (variable_index < this->model.number_variables) { // model variable
         return this->model.variable_lower_bound(variable_index);
      }
      else { // elastic variable in [0, +inf[
         return 0.;
      }
   }

   double l1RelaxedProblem::variable_upper_bound(size_t variable_index) const {
      if (variable_index < this->model.number_variables) { // model variable
         return this->model.variable_upper_bound(variable_index);
      }
      else { // elastic variable in [0, +inf[
         return INF<double>;
      }
   }

   const Vector<size_t>& l1RelaxedProblem::get_fixed_variables() const {
      return this->model.get_fixed_variables();
   }

   double l1RelaxedProblem::constraint_lower_bound(size_t constraint_index) const {
      return this->model.constraint_lower_bound(constraint_index);
   }

   double l1RelaxedProblem::constraint_upper_bound(size_t constraint_index) const {
      return this->model.constraint_upper_bound(constraint_index);
   }

   const Collection<size_t>& l1RelaxedProblem::get_equality_constraints() const {
      return this->model.get_equality_constraints();
   }

   const Collection<size_t>& l1RelaxedProblem::get_inequality_constraints() const {
      return this->model.get_inequality_constraints();
   }

   const Collection<size_t>& l1RelaxedProblem::get_dual_regularization_constraints() const {
      return this->dual_regularization_constraints;
   }

   void l1RelaxedProblem::set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t,
         double)>& elastic_setting_function) const {
      iterate.set_number_variables(this->number_variables);
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
         elastic_setting_function(iterate, constraint_index, elastic_index, -1.);
      }
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
         elastic_setting_function(iterate, constraint_index, elastic_index, 1.);
      }
   }

   void l1RelaxedProblem::set_auxiliary_measure(Iterate& iterate) const {
      iterate.progress.auxiliary = 0.;
      // form the proximal term: zeta/2 ||D_R (x - x_R)||^2
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         double proximal_term = 0.;
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double distance_to_center = iterate.primals[variable_index] - this->proximal_center[variable_index];
            proximal_term += scaling * scaling * distance_to_center * distance_to_center;
         }
         proximal_term *= (this->proximal_coefficient / 2.);
         iterate.progress.auxiliary = proximal_term;
      }
   }

   double l1RelaxedProblem::compute_predicted_auxiliary_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const {
      double predicted_auxiliary_reduction = 0.;
      // form the directional derivative -zeta D_R^2 (x - x_R) and scale it by the step length
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double distance_to_center = current_iterate.primals[variable_index] - this->proximal_center[variable_index];
            predicted_auxiliary_reduction += scaling * scaling * distance_to_center * primal_direction[variable_index];
         }
         predicted_auxiliary_reduction *= step_length * (-this->proximal_coefficient);
      }
      return predicted_auxiliary_reduction;
   }
} // namespace