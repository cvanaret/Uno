// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "OptimizationProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "linear_algebra/MatrixOrder.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Expression.hpp"
#include "tools/Logger.hpp"

namespace uno {
   OptimizationProblem::OptimizationProblem(const Model& model):
         model(model), number_variables(model.number_variables), number_constraints(model.number_constraints),
         primal_regularization_variables(model.number_variables), dual_regularization_constraints(model.number_constraints) {
   }

   OptimizationProblem::OptimizationProblem(const Model& model, size_t number_variables, size_t number_constraints):
         model(model), number_variables(number_variables), number_constraints(number_constraints),
         primal_regularization_variables(model.number_variables), dual_regularization_constraints(model.number_constraints) {
   }

   double OptimizationProblem::get_objective_multiplier() const {
      return 1.;
   }

   void OptimizationProblem::evaluate_constraints(Iterate& iterate, Vector<double>& constraints) const {
      iterate.evaluate_constraints(this->model);
      constraints = iterate.evaluations.constraints;
   }

   void OptimizationProblem::evaluate_objective_gradient(Iterate& iterate, double* objective_gradient) const {
      iterate.evaluate_objective_gradient(this->model);
      for (size_t index: Range(this->number_variables)) {
         objective_gradient[index] = iterate.evaluations.objective_gradient[index];
      }
   }

   size_t OptimizationProblem::number_jacobian_nonzeros() const {
      return this->model.number_jacobian_nonzeros();
   }

   bool OptimizationProblem::has_curvature(const HessianModel& hessian_model) const {
      return hessian_model.has_curvature(this->model);
   }

   size_t OptimizationProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      return hessian_model.number_nonzeros(this->model);
   }

   void OptimizationProblem::compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
         MatrixOrder matrix_order) const {
      this->model.compute_constraint_jacobian_sparsity(row_indices, column_indices, solver_indexing, matrix_order);
   }

   void OptimizationProblem::compute_hessian_sparsity(const HessianModel& hessian_model, int* row_indices,
         int* column_indices, int solver_indexing) const {
      hessian_model.compute_sparsity(this->model, row_indices, column_indices, solver_indexing);
   }

   void OptimizationProblem::evaluate_constraint_jacobian(Iterate& iterate, double* jacobian_values) const {
      this->model.evaluate_constraint_jacobian(iterate.primals, jacobian_values);
   }

   // Lagrangian gradient ∇f(x_k) - ∇c(x_k) y_k - z_k
   // split in two parts: objective contribution and constraints' contribution
   void OptimizationProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient,
         const InequalityHandlingMethod& inequality_handling_method, Iterate& iterate) const {
      lagrangian_gradient.objective_contribution.fill(0.);
      lagrangian_gradient.constraints_contribution.fill(0.);

      // ∇f(x_k)
      this->evaluate_objective_gradient(iterate, lagrangian_gradient.objective_contribution.data());

      // ∇c(x_k) λ_k
      inequality_handling_method.compute_constraint_jacobian_transposed_vector_product(iterate.multipliers.constraints,
         lagrangian_gradient.constraints_contribution);
      lagrangian_gradient.constraints_contribution = -lagrangian_gradient.constraints_contribution;

      // z_k
      for (size_t variable_index: Range(this->number_variables)) {
         lagrangian_gradient.constraints_contribution[variable_index] -= (iterate.multipliers.lower_bounds[variable_index] +
            iterate.multipliers.upper_bounds[variable_index]);
      }
   }

   void OptimizationProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model,
         const Vector<double>& primal_variables, const Multipliers& multipliers, double* hessian_values) const {
      hessian_model.evaluate_hessian(statistics, this->model, primal_variables, this->get_objective_multiplier(),
         multipliers.constraints, hessian_values);
   }

   void OptimizationProblem::compute_hessian_vector_product(HessianModel& hessian_model, const double* x, const double* vector,
         const Multipliers& multipliers, double* result) const {
      hessian_model.compute_hessian_vector_product(this->model, x, vector, this->get_objective_multiplier(),
         multipliers.constraints, result);
   }

   size_t OptimizationProblem::get_number_original_variables() const {
      return this->model.number_variables;
   }

   double OptimizationProblem::variable_lower_bound(size_t variable_index) const {
      return this->model.variable_lower_bound(variable_index);
   }

   double OptimizationProblem::variable_upper_bound(size_t variable_index) const {
      return this->model.variable_upper_bound(variable_index);
   }

   const Vector<size_t>& OptimizationProblem::get_fixed_variables() const {
      return this->model.get_fixed_variables();
   }

   const Collection<size_t>& OptimizationProblem::get_primal_regularization_variables() const {
      return this->primal_regularization_variables;
   }

   double OptimizationProblem::constraint_lower_bound(size_t constraint_index) const {
      return this->model.constraint_lower_bound(constraint_index);
   }

   double OptimizationProblem::constraint_upper_bound(size_t constraint_index) const {
      return this->model.constraint_upper_bound(constraint_index);
   }

   const Collection<size_t>& OptimizationProblem::get_equality_constraints() const {
      return this->model.get_equality_constraints();
   }

   const Collection<size_t>& OptimizationProblem::get_inequality_constraints() const {
      return this->model.get_inequality_constraints();
   }

   const Collection<size_t>& OptimizationProblem::get_dual_regularization_constraints() const {
      return this->dual_regularization_constraints;
   }

   void OptimizationProblem::assemble_primal_dual_direction(const Iterate& /*current_iterate*/, const Vector<double>& /*solution*/,
         Direction& /*direction*/) const {
   }

   double OptimizationProblem::dual_regularization_factor() const {
      return 0.;
   }

   double OptimizationProblem::stationarity_error(const LagrangianGradient<double>& lagrangian_gradient, double objective_multiplier,
         Norm residual_norm) {
      // norm of the scaled Lagrangian gradient
      const auto scaled_lagrangian = objective_multiplier * lagrangian_gradient.objective_contribution + lagrangian_gradient.constraints_contribution;
      return norm(residual_norm, scaled_lagrangian);
   }

   double OptimizationProblem::complementarity_error(const Vector<double>& primals, const Vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const {
      // bound constraints
      const Range variables_range = Range(this->number_variables);
      const VectorExpression variable_complementarity{variables_range, [&](size_t variable_index) {
         assert(variable_index < primals.size());
         assert(variable_index < multipliers.lower_bounds.size());
         assert(variable_index < multipliers.upper_bounds.size());

         if (0. < multipliers.lower_bounds[variable_index]) {
            return multipliers.lower_bounds[variable_index] * (primals[variable_index] - this->variable_lower_bound(variable_index)) - shift_value;
         }
         if (multipliers.upper_bounds[variable_index] < 0.) {
            return multipliers.upper_bounds[variable_index] * (primals[variable_index] - this->variable_upper_bound(variable_index)) - shift_value;
         }
         return 0.;
      }};

      // inequality constraints
      const VectorExpression constraint_complementarity{this->get_inequality_constraints(), [&](size_t constraint_index) {
         assert(constraint_index < constraints.size());
         assert(constraint_index < multipliers.constraints.size());

         if (0. < multipliers.constraints[constraint_index]) { // lower bound
            return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->constraint_lower_bound(constraint_index)) -
               shift_value;
         }
         if (multipliers.constraints[constraint_index] < 0.) { // upper bound
            return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->constraint_upper_bound(constraint_index)) -
               shift_value;
         }
         return 0.;
      }};
      return norm(residual_norm, variable_complementarity, constraint_complementarity);
   }

   SolutionStatus OptimizationProblem::check_first_order_convergence(const Iterate& current_iterate, double primal_tolerance,
         double dual_tolerance) const {
      // evaluate termination conditions based on optimality conditions
      const bool stationarity = (current_iterate.residuals.stationarity / current_iterate.residuals.stationarity_scaling <= dual_tolerance);
      const bool primal_feasibility = (current_iterate.primal_feasibility <= primal_tolerance);
      const bool complementarity = (current_iterate.residuals.complementarity / current_iterate.residuals.complementarity_scaling <= dual_tolerance);

      DEBUG << "\nTermination criteria for primal-dual tolerances = (" << primal_tolerance << ", " << dual_tolerance << "):\n";
      DEBUG << "Stationarity: " << std::boolalpha << stationarity << '\n';
      DEBUG << "Primal feasibility: " << std::boolalpha << primal_feasibility << '\n';
      DEBUG << "Complementarity: " << std::boolalpha << complementarity << '\n';

      if (stationarity && primal_feasibility && 0. < current_iterate.objective_multiplier && complementarity) {
         // feasible regular stationary point
         return SolutionStatus::FEASIBLE_KKT_POINT;
      }
      return SolutionStatus::NOT_OPTIMAL;
   }
} // namespace