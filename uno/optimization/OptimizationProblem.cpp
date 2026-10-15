// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "OptimizationProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Sum.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "tools/Logger.hpp"

namespace uno {
   OptimizationProblem::OptimizationProblem(const Model& model, double residual_scaling_threshold):
         model(model), number_variables(model.number_variables), number_constraints(model.number_constraints),
         residual_scaling_threshold(residual_scaling_threshold), Jv_buffer(this->model.number_constraints) {
   }

   OptimizationProblem::OptimizationProblem(const Model& model, size_t number_variables, size_t number_constraints,
      double residual_scaling_threshold):
         model(model), number_variables(number_variables), number_constraints(number_constraints),
         residual_scaling_threshold(residual_scaling_threshold), Jv_buffer(number_constraints) {
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

   void OptimizationProblem::create_initial_iterate(Iterate& /*iterate*/, Evaluations& /*evaluations*/) const {
      // do nothing
   }

   void OptimizationProblem::postprocess_iterate(Iterate& iterate, Evaluations& /*evaluations*/) const {
      // project the trial iterate onto the bounds to avoid numerical errors
      this->model.project_onto_variable_bounds(iterate.primals);
   }

   void OptimizationProblem::postprocess_solution(Iterate& iterate, Evaluations& evaluations) const {
      this->model.postprocess_solution(iterate, evaluations);
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

   View<const uno_int> OptimizationProblem::get_jacobian_row_indices() const {
      return this->model.get_jacobian_row_indices();
   }

   View<const uno_int> OptimizationProblem::get_jacobian_column_indices() const {
      return this->model.get_jacobian_column_indices();
   }

   void OptimizationProblem::compute_hessian_sparsity(const HessianModel& hessian_model, View<uno_int> row_indices,
         View<uno_int> column_indices, uno_int solver_indexing) const {
      hessian_model.compute_sparsity(row_indices, column_indices, solver_indexing);
   }

   void OptimizationProblem::evaluate_constraints(const Iterate& iterate, View<double> constraints, Evaluations& evaluations) const {
      evaluations.evaluate_constraints(this->model, iterate.primals);
      // copy the constraints at the beginning of the "constraints" vector
      view(constraints, 0, this->model.number_constraints) = evaluations.constraints;
      view(constraints, this->model.number_constraints, constraints.size()).fill(0.);
   }

   void OptimizationProblem::evaluate_objective_gradient(const Iterate& iterate, View<double> objective_gradient,
         Evaluations& evaluations) const {
      evaluations.evaluate_objective_gradient(this->model, iterate.primals);
      // copy the objective gradient at the beginning of the "objective_gradient" vector
      view(objective_gradient, 0, this->model.number_variables) = evaluations.objective_gradient;
      view(objective_gradient, this->model.number_variables, objective_gradient.size()).fill(0.);
   }

   void OptimizationProblem::evaluate_jacobian(const Vector<double>& primals, View<double> jacobian_values, Evaluations& evaluations) const {
      evaluations.evaluate_jacobian(this->model, primals);
      // copy the Jacobian at the beginning of the "jacobian_values" vector
      const size_t number_jacobian_nonzeros = this->model.number_jacobian_nonzeros();
      view(jacobian_values, 0, number_jacobian_nonzeros) = evaluations.jacobian_values;
      view(jacobian_values, number_jacobian_nonzeros, jacobian_values.size()).fill(0.);
   }

   // Lagrangian gradient ∇f(x_k) - ∇c(x_k) y_k - z_k
   void OptimizationProblem::evaluate_lagrangian_gradient(const Iterate& iterate, Evaluations& evaluations,
         Vector<double>& lagrangian_gradient) const {
      this->model.evaluate_lagrangian_gradient(iterate.primals, iterate.multipliers, this->get_objective_multiplier(),
         evaluations, lagrangian_gradient);
   }

   void OptimizationProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model,
         const Vector<double>& primal_variables, const Multipliers& multipliers, View<double> hessian_values) const {
      hessian_model.evaluate_hessian(statistics, primal_variables, this->get_objective_multiplier(),
         multipliers.constraints, hessian_values);
   }

   void OptimizationProblem::add_jacobian_vector_product(View<const double> vector, View<double> result,
         const Evaluations& evaluations) const {
      evaluations.add_jacobian_vector_product(this->model, vector, result);
   }

   void OptimizationProblem::add_jacobian_transposed_vector_product(View<const double> vector, View<double> result,
         const Evaluations& evaluations) const {
      evaluations.add_jacobian_transposed_vector_product(this->model, vector, result);
   }

   void OptimizationProblem::compute_hessian_vector_product(HessianModel& hessian_model, View<const double> x,
         View<const double> vector, const Multipliers& multipliers, View<double> result) const {
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

   Range<> OptimizationProblem::get_primal_regularization_variables() const {
      return Range(this->model.number_variables);
   }

   size_t OptimizationProblem::get_number_bounded_variables() const {
      size_t number_bounded_variables = 0;
      const auto& variables_lower_bounds = this->get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->get_variables_upper_bounds();
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index])) {
            ++number_bounded_variables;
         }
         if (is_finite(variables_upper_bounds[variable_index])) {
            ++number_bounded_variables;
         }
      }
      return number_bounded_variables;
   }

   const std::vector<double>& OptimizationProblem::get_constraints_lower_bounds() const {
      return this->model.get_constraints_lower_bounds();
   }

   const std::vector<double>& OptimizationProblem::get_constraints_upper_bounds() const {
      return this->model.get_constraints_upper_bounds();
   }

   const std::vector<size_t>& OptimizationProblem::get_equality_constraints() const {
      return this->model.get_equality_constraints();
   }

   const std::vector<size_t>& OptimizationProblem::get_inequality_constraints() const {
      return this->model.get_inequality_constraints();
   }

   Inertia OptimizationProblem::get_inertia() const {
      return {this->number_variables, this->number_constraints, 0};
   }

   void OptimizationProblem::assemble_primal_dual_direction(const Iterate& /*current_iterate*/, const Vector<double>& solution,
         Direction& direction) const {
      // form the primal-dual direction
      view(direction.primals, 0, this->number_variables) = view(solution, 0, this->number_variables);
      // retrieve the duals with correct signs (note the sign flip)
      view(direction.multipliers.constraints, 0, this->number_constraints) = -view(solution, this->number_variables,
         this->number_variables + this->number_constraints);
      // set the step lengths (direction unscaled)
      direction.primal_dual_step_length = direction.bound_dual_step_length = 1.;
   }

   void OptimizationProblem::compute_bound_dual_direction(const Vector<double>& /*current_primals*/, const Multipliers& /*current_multipliers*/,
         const Vector<double>& /*direction_primals*/, Multipliers& /*direction_multipliers*/, double& /*bound_dual_step_length*/) const {
      // do nothing
   }

   double OptimizationProblem::dual_regularization_factor() const {
      return 1.;
   }

   double OptimizationProblem::complementarity_error(const Vector<double>& primals, const Vector<double>& constraints,
         const Multipliers& multipliers, Norm residual_norm) const {
      return this->model.complementarity_error(primals, constraints, multipliers, residual_norm);
   }

   SolutionStatus OptimizationProblem::check_first_order_convergence(const Iterate& current_iterate, double primal_tolerance,
         double dual_tolerance) const {
      // evaluate termination conditions based on optimality conditions
      const bool primal_feasibility = (current_iterate.scaled_residuals.primal_infeasibility <= primal_tolerance);
      const bool stationarity = (current_iterate.scaled_residuals.stationarity / current_iterate.scaled_residuals.stationarity_scaling <= dual_tolerance);
      const bool complementarity = (current_iterate.scaled_residuals.complementarity / current_iterate.scaled_residuals.complementarity_scaling <= dual_tolerance);

      DEBUG << "\nTermination criteria for primal-dual tolerances = (" << primal_tolerance << ", " << dual_tolerance << "):\n";
      DEBUG << "Primal infeasibility: " << std::boolalpha << primal_feasibility << '\n';
      DEBUG << "Stationarity: " << std::boolalpha << stationarity << '\n';
      DEBUG << "Complementarity: " << std::boolalpha << complementarity << '\n';

      if (primal_feasibility && stationarity && 0. < current_iterate.objective_multiplier && complementarity) {
         // feasible regular stationary point
         return SolutionStatus::FEASIBLE_KKT_POINT;
      }
      return SolutionStatus::NOT_OPTIMAL;
   }

   double OptimizationProblem::compute_stationarity_scaling(const Multipliers& multipliers) const {
      const size_t total_size = this->get_number_bounded_variables() + this->number_constraints;
      if (total_size == 0) {
         return 1.;
      }
      const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
      const double multiplier_norm = norm_1(multipliers.constraints, multipliers.lower_bounds, multipliers.upper_bounds);
      return std::max(1., multiplier_norm / scaling_factor);
   }

   double OptimizationProblem::compute_complementarity_scaling(const Multipliers& multipliers) const {
      const size_t total_size = this->get_number_bounded_variables();
      if (total_size == 0) {
         return 1.;
      }
      const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
      const double bound_multiplier_norm = norm_1(multipliers.lower_bounds, multipliers.upper_bounds);
      return std::max(1., bound_multiplier_norm / scaling_factor);
   }

   // progress measures

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

   // predicted reductions

   PredictedInfeasibilityReduction OptimizationProblem::build_predicted_infeasibility_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction, Norm norm, Evaluations& current_evaluations) const {
      // α ↦ ‖c(x)‖ − ‖c(x) + ∇c(x)ᵀ(αd)‖
      current_evaluations.evaluate_constraints(this->model, current_iterate.primals);
      current_evaluations.evaluate_jacobian(this->model, current_iterate.primals);
      const double current_constraint_violation = this->model.constraint_violation(current_evaluations.constraints, norm);

      // cache the (expensive) Jacobian–direction product Jd
      this->Jv_buffer.fill(0.);
      current_evaluations.add_jacobian_vector_product(this->model, primal_direction.view(), this->Jv_buffer.view());

      return [this, &current_evaluations, norm, current_constraint_violation](double step_length) {
         return current_constraint_violation - this->model.constraint_violation(current_evaluations.constraints +
            step_length * this->Jv_buffer, norm);
      };
   }

   PredictedObjectiveReduction OptimizationProblem::build_predicted_objective_reduction(
         const Iterate& /*current_iterate*/, const Vector<double>& primal_direction,
         Evaluations& current_evaluations, double hessian_quadratic_form) const {
      // (α, ρ) ↦ −ρ ∇f(x)ᵀ(αd) − α²/2 dᵀHd   (ρ scales only the gradient term; H already carries ρ)
      const double directional_derivative = dot(view(primal_direction, 0, this->model.number_variables),
         current_evaluations.objective_gradient);
      return [directional_derivative, hessian_quadratic_form](double step_length, double objective_multiplier, bool use_curvature_information) {
         double reduction = -objective_multiplier * step_length * directional_derivative;
         if (use_curvature_information) {
            reduction -= step_length * step_length / 2. * hessian_quadratic_form;
         }
         return reduction;
      };
   }

   PredictedAuxiliaryReduction OptimizationProblem::build_predicted_auxiliary_reduction(
         const Iterate& /*current_iterate*/, const Vector<double>& /*primal_direction*/) const {
      return [](double /*step_length*/) { return 0.; };
   }
} // namespace