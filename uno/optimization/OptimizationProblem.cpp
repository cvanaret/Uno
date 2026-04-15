// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "OptimizationProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "linear_algebra/MatrixOrder.hpp"
#include "linear_algebra/VectorView.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/UnaryNegation.hpp"
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

   std::unique_ptr<OptimizationProblem> OptimizationProblem::clone() const {
      return std::make_unique<OptimizationProblem>(*this);
   }

   double OptimizationProblem::get_objective_multiplier() const {
      return 1.;
   }

   bool OptimizationProblem::has_inequality_constraints() const {
      return this->model.has_inequality_constraints();
   }

   bool OptimizationProblem::has_bound_constraints() const {
      return this->model.has_bound_constraints();
   }

   void OptimizationProblem::generate_initial_iterate(Iterate& /*initial_iterate*/, Evaluations& /*evaluations*/) const {
      // do nothing
   }

   void OptimizationProblem::postprocess_iterate(Iterate& /*iterate*/) const {
      // do nothing
   }

   size_t OptimizationProblem::number_jacobian_nonzeros() const {
      return this->model.number_jacobian_nonzeros();
   }

   bool OptimizationProblem::has_curvature(const HessianModel& hessian_model) const {
      return hessian_model.has_curvature();
   }

   size_t OptimizationProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      return hessian_model.number_nonzeros();
   }

   void OptimizationProblem::compute_jacobian_sparsity(uno_int *row_indices, uno_int *column_indices, uno_int row_offset,
         uno_int column_offset, uno_int solver_indexing, MatrixOrder matrix_order) const {
      this->model.compute_jacobian_sparsity(row_indices, column_indices, row_offset, column_offset, solver_indexing, matrix_order);
   }

   void OptimizationProblem::compute_hessian_sparsity(const HessianModel& hessian_model, uno_int *row_indices,
         uno_int *column_indices, uno_int solver_indexing) const {
      hessian_model.compute_sparsity(row_indices, column_indices, solver_indexing);
   }

   void OptimizationProblem::evaluate_constraints(const Iterate& iterate, double* constraints, Evaluations& evaluations) const {
      evaluations.evaluate_constraints(this->model, iterate.primals);
      for (size_t index: Range(this->number_constraints)) {
         constraints[index] = evaluations.constraints[index];
      }
   }

   // warning: adds to objective_gradient (objective_gradient must be reset prior, if necessary)
   void OptimizationProblem::evaluate_objective_gradient(const Iterate& iterate, double* objective_gradient, Evaluations& evaluations) const {
      evaluations.evaluate_objective_gradient(this->model, iterate.primals);
      for (size_t index: Range(this->number_variables)) {
         objective_gradient[index] += evaluations.objective_gradient[index];
      }
   }

   void OptimizationProblem::evaluate_jacobian(const Vector<double>& primals, double* jacobian_values, Evaluations& evaluations) const {
      evaluations.evaluate_jacobian(this->model, primals);
      for (size_t nonzero_index: Range(this->model.number_jacobian_nonzeros())) {
         jacobian_values[nonzero_index] = evaluations.jacobian_values[nonzero_index];
      }
   }

   // Lagrangian gradient ∇f(x_k) - ∇c(x_k) y_k - z_k
   void OptimizationProblem::evaluate_lagrangian_gradient(const Iterate& iterate, Evaluations& evaluations,
         Vector<double>& lagrangian_gradient) const {
      this->model.evaluate_lagrangian_gradient(iterate.primals, iterate.multipliers, 1., evaluations, lagrangian_gradient);
   }

   void OptimizationProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model,
         const Vector<double>& primal_variables, const Multipliers& multipliers, double* hessian_values) const {
      hessian_model.evaluate_hessian(statistics, primal_variables, this->get_objective_multiplier(),
         multipliers.constraints, hessian_values);
   }

   void OptimizationProblem::compute_jacobian_vector_product(const double* vector, double* result,
         const Evaluations& evaluations) const {
      evaluations.compute_jacobian_vector_product(this->model, vector, result);
   }

   void OptimizationProblem::compute_jacobian_transposed_vector_product(const double* vector, double* result,
         const Evaluations& evaluations) const {
      evaluations.compute_jacobian_transposed_vector_product(this->model, vector, result);
   }

   void OptimizationProblem::compute_hessian_vector_product(HessianModel& hessian_model, const double* x, const double* vector,
         const Multipliers& multipliers, double* result) const {
      hessian_model.compute_hessian_vector_product(x, vector, this->get_objective_multiplier(), multipliers.constraints, result);
   }

   size_t OptimizationProblem::get_number_original_variables() const {
      return this->model.number_variables;
   }

   const std::vector<double>& OptimizationProblem::get_variables_lower_bounds() const {
      return this->model.get_variables_lower_bounds();
   }

   const std::vector<double>& OptimizationProblem::get_variables_upper_bounds() const {
      return this->model.get_variables_upper_bounds();
   }

   const Vector<size_t>& OptimizationProblem::get_fixed_variables() const {
      return this->model.get_fixed_variables();
   }

   const Collection<size_t>& OptimizationProblem::get_primal_regularization_variables() const {
      return this->primal_regularization_variables;
   }

   const std::vector<double>& OptimizationProblem::get_constraints_lower_bounds() const {
      return this->model.get_constraints_lower_bounds();
   }

   const std::vector<double>& OptimizationProblem::get_constraints_upper_bounds() const {
      return this->model.get_constraints_upper_bounds();
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

   Inertia OptimizationProblem::get_inertia() const {
      return {this->number_variables, this->number_constraints, 0};
   }

   void OptimizationProblem::assemble_primal_dual_direction(const Iterate& /*current_iterate*/, const Vector<double>& solution,
         Direction& direction) const {
      // form the primal-dual direction
      view(direction.primals, 0, this->number_variables) = view(solution, 0, this->number_variables);
      // retrieve the duals with correct signs (note the sign flip)
      direction.multipliers.constraints = -view(solution, this->number_variables, this->number_variables + this->number_constraints);
   }

   double OptimizationProblem::dual_regularization_factor() const {
      return 0.;
   }
   
   double OptimizationProblem::complementarity_error(const Vector<double>& primals, const Vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const {
      // bound constraints
      const Range variables_range = Range(this->number_variables);
      const auto& variables_lower_bounds = this->get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->get_variables_upper_bounds();
      const VectorExpression variable_complementarity{variables_range, [&](size_t variable_index) {
         assert(variable_index < primals.size());
         assert(variable_index < multipliers.lower_bounds.size());
         assert(variable_index < multipliers.upper_bounds.size());

         if (0. < multipliers.lower_bounds[variable_index]) {
            return multipliers.lower_bounds[variable_index] * (primals[variable_index] - variables_lower_bounds[variable_index]) - shift_value;
         }
         if (multipliers.upper_bounds[variable_index] < 0.) {
            return multipliers.upper_bounds[variable_index] * (primals[variable_index] - variables_upper_bounds[variable_index]) - shift_value;
         }
         return 0.;
      }};

      // inequality constraints
      const auto& constraints_lower_bounds = this->model.get_constraints_lower_bounds();
      const auto& constraints_upper_bounds = this->model.get_constraints_upper_bounds();
      const VectorExpression constraint_complementarity{this->get_inequality_constraints(), [&](size_t constraint_index) {
         assert(constraint_index < constraints.size());
         assert(constraint_index < multipliers.constraints.size());

         if (0. < multipliers.constraints[constraint_index]) { // lower bound
            return multipliers.constraints[constraint_index] * (constraints[constraint_index] - constraints_lower_bounds[constraint_index]) -
               shift_value;
         }
         if (multipliers.constraints[constraint_index] < 0.) { // upper bound
            return multipliers.constraints[constraint_index] * (constraints[constraint_index] - constraints_upper_bounds[constraint_index]) -
               shift_value;
         }
         return 0.;
      }};
      return norm(residual_norm, variable_complementarity, constraint_complementarity);
   }

   double OptimizationProblem::compute_centrality_error(const Vector<double>& primals, const Multipliers& multipliers,
         double shift) const {
      const Range variables_range = Range(this->number_variables);
      const auto& variables_lower_bounds = this->get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->get_variables_upper_bounds();
      const VectorExpression shifted_bound_complementarity{variables_range, [&](size_t variable_index) {
         double result = 0.;
         if (0. < multipliers.lower_bounds[variable_index]) { // lower bound
            result = std::max(result, std::abs(multipliers.lower_bounds[variable_index] *
               (primals[variable_index] - variables_lower_bounds[variable_index]) - shift));
         }
         if (multipliers.upper_bounds[variable_index] < 0.) { // upper bound
            result = std::max(result, std::abs(multipliers.upper_bounds[variable_index] *
               (primals[variable_index] - variables_upper_bounds[variable_index]) - shift));
         }
         return result;
      }};
      return norm_inf(shifted_bound_complementarity); // TODO use a generic norm
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

   // infeasibility measure: constraint violation
   void OptimizationProblem::set_infeasibility_measure(Iterate& iterate, Evaluations& evaluations, Norm norm) const {
      evaluations.evaluate_constraints(this->model, iterate.primals);
      iterate.progress.infeasibility = this->model.constraint_violation(evaluations.constraints, norm);
   }

   // objective measure: scaled objective
   void OptimizationProblem::set_objective_measure(Iterate& iterate, Evaluations& evaluations) const {
      evaluations.evaluate_objective(this->model, iterate.primals);
      const double objective = evaluations.objective;
      iterate.progress.objective = [=](double objective_multiplier) {
         return objective_multiplier * objective;
      };
   }

   void OptimizationProblem::set_auxiliary_measure(Iterate& iterate) const {
      iterate.progress.auxiliary = 0.;
   }

   double OptimizationProblem::compute_predicted_auxiliary_reduction(const Iterate& /*current_iterate*/,
         const Vector<double>& /*primal_direction*/, double /*step_length*/) const {
      return 0.;
   }
} // namespace