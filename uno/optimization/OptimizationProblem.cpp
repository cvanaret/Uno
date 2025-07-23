// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "OptimizationProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Expression.hpp"

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

   void OptimizationProblem::evaluate_objective_gradient(Iterate& iterate, Vector<double>& objective_gradient) const {
      iterate.evaluate_objective_gradient(this->model);
      // TODO change this
      objective_gradient = iterate.evaluations.objective_gradient;
   }

   void OptimizationProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
      iterate.evaluate_constraints(this->model);
      constraints = iterate.evaluations.constraints;
   }

   void OptimizationProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
      iterate.evaluate_constraint_jacobian(this->model);
      // TODO change this
      constraint_jacobian = iterate.evaluations.constraint_jacobian;
   }

   void OptimizationProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, SymmetricMatrix<size_t, double>& hessian) const {
      hessian_model.evaluate_hessian(statistics, this->model, primal_variables, this->get_objective_multiplier(), multipliers.constraints, hessian);
   }

   void OptimizationProblem::compute_hessian_vector_product(HessianModel& hessian_model, const double* vector, const Multipliers& multipliers,
         double* result) const {
      hessian_model.compute_hessian_vector_product(this->model, vector, this->get_objective_multiplier(), multipliers.constraints, result);
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

   const Collection<size_t>& OptimizationProblem::get_lower_bounded_variables() const {
      return this->model.get_lower_bounded_variables();
   }

   const Collection<size_t>& OptimizationProblem::get_upper_bounded_variables() const {
      return this->model.get_upper_bounded_variables();
   }

   const Collection<size_t>& OptimizationProblem::get_single_lower_bounded_variables() const {
      return this->model.get_single_lower_bounded_variables();
   }

   const Collection<size_t>& OptimizationProblem::get_single_upper_bounded_variables() const {
      return this->model.get_single_upper_bounded_variables();
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

   size_t OptimizationProblem::number_jacobian_nonzeros() const {
      return this->model.number_jacobian_nonzeros();
   }

   size_t OptimizationProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      return hessian_model.number_nonzeros(this->model);
   }

   void OptimizationProblem::assemble_primal_dual_direction(const Iterate& /*current_iterate*/, const Multipliers& /*current_multipliers*/,
         const Vector<double>& /*solution*/, Direction& /*direction*/) const {
      // do nothing
   }

   double OptimizationProblem::stationarity_error(const LagrangianGradient<double>& lagrangian_gradient, double objective_multiplier,
         Norm residual_norm) {
      // norm of the scaled Lagrangian gradient
      const auto scaled_lagrangian = objective_multiplier * lagrangian_gradient.objective_contribution + lagrangian_gradient.constraints_contribution;
      return norm(residual_norm, scaled_lagrangian);
   }

   // Lagrangian gradient split in two parts: objective contribution and constraints' contribution
   void OptimizationProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate,
         const Multipliers& multipliers) const {
      lagrangian_gradient.objective_contribution.fill(0.);
      lagrangian_gradient.constraints_contribution.fill(0.);

      // objective gradient
      lagrangian_gradient.objective_contribution = iterate.evaluations.objective_gradient;

      // constraints
      for (size_t constraint_index: Range(this->number_constraints)) {
         if (multipliers.constraints[constraint_index] != 0.) {
            for (auto [variable_index, derivative]: iterate.evaluations.constraint_jacobian[constraint_index]) {
               lagrangian_gradient.constraints_contribution[variable_index] -= multipliers.constraints[constraint_index] * derivative;
            }
         }
      }

      // bound constraints of original variables
      for (size_t variable_index: Range(this->number_variables)) {
         lagrangian_gradient.constraints_contribution[variable_index] -= (multipliers.lower_bounds[variable_index] +
                                                                          multipliers.upper_bounds[variable_index]);
      }
   }

   double OptimizationProblem::complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const {
      // bound constraints
      const Range variables_range = Range(std::min(this->number_variables, primals.size()));
      const VectorExpression variable_complementarity{variables_range, [&](size_t variable_index) {
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
         if (0. < multipliers.constraints[constraint_index]) { // lower bound
            return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->constraint_lower_bound(constraint_index)) -
               shift_value;
         }
         else if (multipliers.constraints[constraint_index] < 0.) { // upper bound
            return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->constraint_upper_bound(constraint_index)) -
               shift_value;
         }
         return 0.;
      }};
      return norm(residual_norm, variable_complementarity, constraint_complementarity);
   }

   double OptimizationProblem::dual_regularization_factor() const {
      return 0.;
   }
} // namespace