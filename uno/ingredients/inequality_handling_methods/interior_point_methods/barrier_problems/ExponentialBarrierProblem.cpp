// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ExponentialBarrierProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   ExponentialBarrierProblem::ExponentialBarrierProblem(const OptimizationProblem& problem, double barrier_parameter,
      const InteriorPointParameters& parameters):
         // no slacks: as many constraints as the number of equality constraints of the problem
         BarrierProblem(problem.model, problem.number_variables, problem.get_equality_constraints().size()),
         reformulated_problem(problem), barrier_parameter(barrier_parameter), parameters(parameters) {
      if (!this->reformulated_problem.get_inequality_constraints().empty()) {
         throw std::runtime_error("ExponentialBarrierProblem does not support inequality constraints yet");
      }
   }

   double ExponentialBarrierProblem::get_objective_multiplier() const {
      return this->reformulated_problem.get_objective_multiplier();
   }

   void ExponentialBarrierProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
      this->reformulated_problem.evaluate_constraints(iterate, constraints);
   }

   void ExponentialBarrierProblem::evaluate_objective_gradient(Iterate& iterate, double* objective_gradient) const {
      this->reformulated_problem.evaluate_objective_gradient(iterate, objective_gradient);
   }

   size_t ExponentialBarrierProblem::number_jacobian_nonzeros() const {
      return this->reformulated_problem.number_jacobian_nonzeros();
   }

   bool ExponentialBarrierProblem::has_curvature(const HessianModel& hessian_model) const {
      return hessian_model.has_curvature(this->model);
   }

   size_t ExponentialBarrierProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      size_t number_nonzeros = this->reformulated_problem.number_hessian_nonzeros(hessian_model);
      return number_nonzeros;
   }

   void ExponentialBarrierProblem::compute_constraint_jacobian_sparsity(size_t* row_indices, size_t* column_indices,
         size_t solver_indexing, MatrixOrder matrix_order) const {
      this->reformulated_problem.compute_constraint_jacobian_sparsity(row_indices, column_indices, solver_indexing, matrix_order);
   }

   void ExponentialBarrierProblem::compute_hessian_sparsity(const HessianModel& hessian_model, size_t* row_indices,
         size_t* column_indices, size_t solver_indexing) const {
      // original Lagrangian Hessian
      this->reformulated_problem.compute_hessian_sparsity(hessian_model, row_indices, column_indices, solver_indexing);
   }

   void ExponentialBarrierProblem::evaluate_constraint_jacobian(Iterate& iterate, double* jacobian_values) const {
      this->reformulated_problem.evaluate_constraint_jacobian(iterate, jacobian_values);
   }

   void ExponentialBarrierProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient,
         const InequalityHandlingMethod& inequality_handling_method, Iterate& iterate) const {
      this->reformulated_problem.evaluate_lagrangian_gradient(lagrangian_gradient, inequality_handling_method, iterate);
   }
   
   void ExponentialBarrierProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model,
         const Vector<double>& primal_variables, const Multipliers& multipliers, double* hessian_values) const {
      // original Lagrangian Hessian
      this->reformulated_problem.evaluate_lagrangian_hessian(statistics, hessian_model, primal_variables, multipliers, hessian_values);
   }

   void ExponentialBarrierProblem::compute_hessian_vector_product(HessianModel& hessian_model, const double* vector,
         const Multipliers& multipliers, double* result) const {
      // original Lagrangian Hessian
      this->reformulated_problem.compute_hessian_vector_product(hessian_model, vector, multipliers, result);
   }

   double ExponentialBarrierProblem::variable_lower_bound(size_t /*variable_index*/) const {
      return -INF<double>;
   }

   double ExponentialBarrierProblem::variable_upper_bound(size_t /*variable_index*/) const {
      return INF<double>;
   }

   const Collection<size_t>& ExponentialBarrierProblem::get_lower_bounded_variables() const {
      return this->reformulated_problem.get_lower_bounded_variables();
   }

   const Collection<size_t>& ExponentialBarrierProblem::get_upper_bounded_variables() const {
      return this->reformulated_problem.get_upper_bounded_variables();
   }

   const Collection<size_t>& ExponentialBarrierProblem::get_single_lower_bounded_variables() const {
      return this->reformulated_problem.get_single_lower_bounded_variables();
   }

   const Collection<size_t>& ExponentialBarrierProblem::get_single_upper_bounded_variables() const {
      return this->reformulated_problem.get_single_upper_bounded_variables();
   }

   const Vector<size_t>& ExponentialBarrierProblem::get_fixed_variables() const {
      return this->fixed_variables;
   }

   const Collection<size_t>& ExponentialBarrierProblem::get_primal_regularization_variables() const {
      return this->reformulated_problem.get_primal_regularization_variables();
   }

   double ExponentialBarrierProblem::constraint_lower_bound(size_t constraint_index) const {
      return this->reformulated_problem.constraint_lower_bound(constraint_index);
   }

   double ExponentialBarrierProblem::constraint_upper_bound(size_t constraint_index) const {
      return this->reformulated_problem.constraint_upper_bound(constraint_index);
   }

   const Collection<size_t>& ExponentialBarrierProblem::get_equality_constraints() const {
      return this->reformulated_problem.get_equality_constraints();
   }

   const Collection<size_t>& ExponentialBarrierProblem::get_inequality_constraints() const {
      return this->inequality_constraints;
   }

   const Collection<size_t>& ExponentialBarrierProblem::get_dual_regularization_constraints() const {
      if (this->reformulated_problem.get_dual_regularization_constraints().empty()) {
         // this is an indication that the constraints (if there is any) were already regularized in a previous
         // reformulation (e.g. l1 relaxation). In that case, we stick to an empty set
         return this->reformulated_problem.get_dual_regularization_constraints();
      }
      // otherwise, we pick the set of equality constraints
      return this->reformulated_problem.get_equality_constraints();
   }

   void ExponentialBarrierProblem::assemble_primal_dual_direction(const Iterate& current_iterate, const Vector<double>& solution,
         Direction& direction) const {
      // form the primal-dual direction
      direction.primals = view(solution, 0, this->reformulated_problem.number_variables);
      // retrieve the duals with correct signs (note the minus sign)
      direction.multipliers.constraints = view(-solution, this->reformulated_problem.number_variables,
         this->reformulated_problem.number_variables + this->reformulated_problem.number_constraints);

      // "fraction-to-boundary" rule for primal variables and constraints multipliers
      const double tau = std::max(this->parameters.tau_min, 1. - this->barrier_parameter);
      const double primal_step_length = ExponentialBarrierProblem::primal_fraction_to_boundary(current_iterate.primals,
         direction.primals, tau);
      const double bound_dual_step_length = ExponentialBarrierProblem::dual_fraction_to_boundary(current_iterate.multipliers,
         direction.multipliers, tau);
      DEBUG << "Fraction-to-boundary rules:\n";
      DEBUG << "primal step length = " << primal_step_length << '\n';
      DEBUG << "bound dual step length = " << bound_dual_step_length << "\n\n";
      // scale the primal-dual variables
      direction.primals.scale(primal_step_length);
      direction.multipliers.constraints.scale(primal_step_length);
      direction.multipliers.lower_bounds.scale(bound_dual_step_length);
      direction.multipliers.upper_bounds.scale(bound_dual_step_length);
   }

   double ExponentialBarrierProblem::push_variable_to_interior(double variable_value, double /*lower_bound*/, double /*upper_bound*/) const {
      return variable_value;
   }

   void ExponentialBarrierProblem::set_auxiliary_measure(Iterate& iterate) const {
      // auxiliary measure: barrier terms
      double barrier_terms = 0.;
      iterate.progress.auxiliary = barrier_terms;
   }

   double ExponentialBarrierProblem::complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const {
      return this->reformulated_problem.complementarity_error(primals, constraints, multipliers, shift_value, residual_norm);
   }

   double ExponentialBarrierProblem::dual_regularization_factor() const {
      return std::pow(this->barrier_parameter, this->parameters.dual_regularization_exponent);
   }

   // protected member functions

   double ExponentialBarrierProblem::compute_barrier_term_directional_derivative(const Iterate& /*current_iterate*/,
         const Vector<double>& /*primal_direction*/) const {
      double directional_derivative = 0.;
      return directional_derivative;
   }

   // TODO use a single function for primal and dual fraction-to-boundary rules
   double ExponentialBarrierProblem::primal_fraction_to_boundary(const Vector<double>& current_primals,
         const Vector<double>& primal_direction, double tau) const {
      double step_length = 1.;
      for (const size_t variable_index: this->reformulated_problem.get_lower_bounded_variables()) {
         if (primal_direction[variable_index] < 0.) {
            const double distance = -tau * (current_primals[variable_index] - this->reformulated_problem.variable_lower_bound(variable_index)) / primal_direction[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      for (const size_t variable_index: this->reformulated_problem.get_upper_bounded_variables()) {
         if (0. < primal_direction[variable_index]) {
            const double distance = -tau * (current_primals[variable_index] - this->reformulated_problem.variable_upper_bound(variable_index)) / primal_direction[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      assert(0. < step_length && step_length <= 1. && "The primal fraction-to-boundary step length is not in (0, 1]");
      return step_length;
   }

   double ExponentialBarrierProblem::dual_fraction_to_boundary(const Multipliers& current_multipliers,
         const Multipliers& direction_multipliers, double tau) const {
      double step_length = 1.;
      for (const size_t variable_index: this->reformulated_problem.get_lower_bounded_variables()) {
         if (direction_multipliers.lower_bounds[variable_index] < 0.) {
            const double distance = -tau * current_multipliers.lower_bounds[variable_index] / direction_multipliers.lower_bounds[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      for (const size_t variable_index: this->reformulated_problem.get_upper_bounded_variables()) {
         if (0. < direction_multipliers.upper_bounds[variable_index]) {
            const double distance = -tau * current_multipliers.upper_bounds[variable_index] / direction_multipliers.upper_bounds[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      assert(0. < step_length && step_length <= 1. && "The dual fraction-to-boundary step length is not in (0, 1]");
      return step_length;
   }

   void ExponentialBarrierProblem::postprocess_iterate(Iterate& /*iterate*/) const {
      // do nothing
   }

   double ExponentialBarrierProblem::compute_centrality_error(const Vector<double>& primals,
         const Multipliers& multipliers, double barrier_parameter) const {
      const Range variables_range = Range(this->reformulated_problem.number_variables);
      const VectorExpression shifted_bound_complementarity{variables_range, [&](size_t variable_index) {
         double result = 0.;
         if (0. < multipliers.lower_bounds[variable_index]) { // lower bound
            result = std::max(result, std::abs(multipliers.lower_bounds[variable_index] *
               (primals[variable_index] - this->reformulated_problem.variable_lower_bound(variable_index)) - barrier_parameter));
         }
         if (multipliers.upper_bounds[variable_index] < 0.) { // upper bound
            result = std::max(result, std::abs(multipliers.upper_bounds[variable_index] *
               (primals[variable_index] - this->reformulated_problem.variable_upper_bound(variable_index)) - barrier_parameter));
         }
         return result;
      }};
      return norm_inf(shifted_bound_complementarity); // TODO use a generic norm
   }
} // namespace