// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalDualInteriorPointProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   PrimalDualInteriorPointProblem::PrimalDualInteriorPointProblem(const OptimizationProblem& problem, double barrier_parameter,
      const InteriorPointParameters &parameters):
         OptimizationProblem(problem.model, problem.number_variables, problem.number_constraints),
         first_reformulation(problem), barrier_parameter(barrier_parameter),
         parameters(parameters), equality_constraints(problem.number_constraints) { }

   double PrimalDualInteriorPointProblem::get_objective_multiplier() const {
      return this->first_reformulation.get_objective_multiplier();
   }

   void PrimalDualInteriorPointProblem::evaluate_constraints(Iterate& iterate, Vector<double>& constraints) const {
      this->first_reformulation.evaluate_constraints(iterate, constraints);
   }

   void PrimalDualInteriorPointProblem::evaluate_objective_gradient(Iterate& iterate, double* objective_gradient) const {
      this->first_reformulation.evaluate_objective_gradient(iterate, objective_gradient);

      // barrier terms
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         double barrier_term = 0.;
         if (is_finite(this->first_reformulation.variable_lower_bound(variable_index))) { // lower bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - this->first_reformulation.variable_lower_bound(variable_index));
            // damping
            if (!is_finite(this->first_reformulation.variable_upper_bound(variable_index))) {
               barrier_term += this->parameters.damping_factor * this->barrier_parameter;
            }
         }
         if (is_finite(this->first_reformulation.variable_upper_bound(variable_index))) { // upper bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - this->first_reformulation.variable_upper_bound(variable_index));
            // damping
            if (!is_finite(this->first_reformulation.variable_lower_bound(variable_index))) {
               barrier_term -= this->parameters.damping_factor * this->barrier_parameter;
            }
         }
         objective_gradient[variable_index] += barrier_term;
      }
   }

   void PrimalDualInteriorPointProblem::compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices,
         int solver_indexing, MatrixOrder matrix_order) const {
      this->first_reformulation.compute_constraint_jacobian_sparsity(row_indices, column_indices, solver_indexing, matrix_order);
   }

   void PrimalDualInteriorPointProblem::compute_hessian_sparsity(const HessianModel& hessian_model, int* row_indices,
         int* column_indices, int solver_indexing) const {
      // original Lagrangian Hessian
      this->first_reformulation.compute_hessian_sparsity(hessian_model, row_indices, column_indices, solver_indexing);

      // diagonal barrier terms
      size_t current_index = this->first_reformulation.number_hessian_nonzeros(hessian_model);
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const bool finite_lower_bound = is_finite(this->first_reformulation.variable_lower_bound(variable_index));
         const bool finite_upper_bound = is_finite(this->first_reformulation.variable_upper_bound(variable_index));
         if (finite_lower_bound || finite_upper_bound) {
            row_indices[current_index] = static_cast<int>(variable_index) + solver_indexing;
            column_indices[current_index] = static_cast<int>(variable_index) + solver_indexing;
            ++current_index;
         }
      }
   }

   size_t PrimalDualInteriorPointProblem::number_jacobian_nonzeros() const {
      return this->first_reformulation.number_jacobian_nonzeros();
   }

   bool PrimalDualInteriorPointProblem::has_curvature(const HessianModel& hessian_model) const {
      if (hessian_model.has_curvature(this->model)) {
         return true;
      }
      else {
         // barrier terms
         for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
            if (is_finite(this->first_reformulation.variable_lower_bound(variable_index)) || is_finite(this->first_reformulation.variable_upper_bound(variable_index))) {
               return true;
            }
         }
         return false;
      }
   }

   size_t PrimalDualInteriorPointProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      size_t number_nonzeros = this->first_reformulation.number_hessian_nonzeros(hessian_model);
      // barrier contribution: original variables
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const double lower_bound = this->first_reformulation.variable_lower_bound(variable_index);
         const double upper_bound = this->first_reformulation.variable_upper_bound(variable_index);
         if (is_finite(lower_bound) || is_finite(upper_bound)) {
            ++number_nonzeros;
         }
      }
      return number_nonzeros;
   }

   void PrimalDualInteriorPointProblem::evaluate_constraint_jacobian(Iterate& iterate, double* jacobian_values) const {
      this->first_reformulation.evaluate_constraint_jacobian(iterate, jacobian_values);
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient,
         const InequalityHandlingMethod& inequality_handling_method, Iterate& iterate) const {
      this->first_reformulation.evaluate_lagrangian_gradient(lagrangian_gradient, inequality_handling_method, iterate);

      // barrier terms
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         double barrier_term = 0.;
         if (is_finite(this->first_reformulation.variable_lower_bound(variable_index))) { // lower bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - this->first_reformulation.variable_lower_bound(variable_index));
            // damping
            if (!is_finite(this->first_reformulation.variable_upper_bound(variable_index))) {
               barrier_term += this->parameters.damping_factor * this->barrier_parameter;
            }
         }
         if (is_finite(this->first_reformulation.variable_upper_bound(variable_index))) { // upper bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - this->first_reformulation.variable_upper_bound(variable_index));
            // damping
            if (!is_finite(this->first_reformulation.variable_lower_bound(variable_index))) {
               barrier_term -= this->parameters.damping_factor * this->barrier_parameter;
            }
         }
         // the objective contribution of the Lagrangian gradient may be scaled. Barrier terms go into the constraint contribution
         lagrangian_gradient.constraints_contribution[variable_index] += barrier_term;
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, double* hessian_values) const {
      // original Lagrangian Hessian
      this->first_reformulation.evaluate_lagrangian_hessian(statistics, hessian_model, primal_variables, multipliers, hessian_values);

      // barrier terms
      size_t nonzero_index = this->first_reformulation.number_hessian_nonzeros(hessian_model);
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const bool finite_lower_bound = is_finite(this->first_reformulation.variable_lower_bound(variable_index));
         const bool finite_upper_bound = is_finite(this->first_reformulation.variable_upper_bound(variable_index));
         if (finite_lower_bound || finite_upper_bound) {
            double diagonal_barrier_term = 0.;
            if (finite_lower_bound) {
               const double distance_to_bound = primal_variables[variable_index] - this->first_reformulation.variable_lower_bound(variable_index);
               diagonal_barrier_term += multipliers.lower_bounds[variable_index] / distance_to_bound;
            }
            if (finite_upper_bound) {
               const double distance_to_bound = primal_variables[variable_index] - this->first_reformulation.variable_upper_bound(variable_index);
               diagonal_barrier_term += multipliers.upper_bounds[variable_index] / distance_to_bound;
            }
            hessian_values[nonzero_index] = diagonal_barrier_term;
            ++nonzero_index;
         }
      }
   }

   void PrimalDualInteriorPointProblem::compute_hessian_vector_product(HessianModel& hessian_model, const double* x,
         const double* vector, const Multipliers& multipliers, double* result) const {
      // original Lagrangian Hessian
      this->first_reformulation.compute_hessian_vector_product(hessian_model, x, vector, multipliers, result);

      // barrier terms
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const bool finite_lower_bound = is_finite(this->first_reformulation.variable_lower_bound(variable_index));
         const bool finite_upper_bound = is_finite(this->first_reformulation.variable_upper_bound(variable_index));
         if (finite_lower_bound || finite_upper_bound) {
            double diagonal_barrier_term = 0.;
            if (finite_lower_bound) { // lower bounded
               const double distance_to_bound = vector[variable_index] - this->first_reformulation.variable_lower_bound(variable_index);
               diagonal_barrier_term += multipliers.lower_bounds[variable_index] / distance_to_bound;
            }
            if (finite_upper_bound) { // upper bounded
               const double distance_to_bound = vector[variable_index] - this->first_reformulation.variable_upper_bound(variable_index);
               diagonal_barrier_term += multipliers.upper_bounds[variable_index] / distance_to_bound;
            }
            result[variable_index] += diagonal_barrier_term * vector[variable_index];
         }
      }
   }

   double PrimalDualInteriorPointProblem::variable_lower_bound(size_t /*variable_index*/) const {
      return -INF<double>;
   }

   double PrimalDualInteriorPointProblem::variable_upper_bound(size_t /*variable_index*/) const {
      return INF<double>;
   }

   const Vector<size_t>& PrimalDualInteriorPointProblem::get_fixed_variables() const {
      return this->fixed_variables;
   }

   double PrimalDualInteriorPointProblem::constraint_lower_bound(size_t /*constraint_index*/) const {
      return 0.;
   }

   double PrimalDualInteriorPointProblem::constraint_upper_bound(size_t /*constraint_index*/) const {
      return 0.;
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_equality_constraints() const {
      return this->equality_constraints;
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_inequality_constraints() const {
      return this->inequality_constraints;
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_dual_regularization_constraints() const {
      if (this->first_reformulation.get_dual_regularization_constraints().empty()) {
         // this is an indication that the constraints (if there is any) were already regularized in a previous
         // reformulation (e.g. l1 relaxation). In that case, we stick to an empty set
         return this->first_reformulation.get_dual_regularization_constraints();
      }
      // otherwise, we pick the set of equality constraints, since the inequality constraints have slacks
      return this->first_reformulation.get_equality_constraints();
   }

   void PrimalDualInteriorPointProblem::assemble_primal_dual_direction(const Iterate& current_iterate, const Vector<double>& solution,
         Direction& direction) const {
      // form the primal-dual direction
      direction.primals = view(solution, 0, this->first_reformulation.number_variables);
      // retrieve the duals with correct signs (note the minus sign)
      direction.multipliers.constraints = view(-solution, this->first_reformulation.number_variables,
         this->first_reformulation.number_variables + this->first_reformulation.number_constraints);
      this->compute_bound_dual_direction(current_iterate, direction);

      // "fraction-to-boundary" rule for primal variables and constraints multipliers
      const double tau = std::max(this->parameters.tau_min, 1. - this->barrier_parameter);
      const double primal_step_length = PrimalDualInteriorPointProblem::primal_fraction_to_boundary(current_iterate.primals,
         direction.primals, tau);
      const double bound_dual_step_length = PrimalDualInteriorPointProblem::dual_fraction_to_boundary(current_iterate.multipliers,
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

   void PrimalDualInteriorPointProblem::set_auxiliary_measure(Iterate& iterate) const {
      // auxiliary measure: barrier terms
      double barrier_terms = 0.;
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const double lower_bound = this->first_reformulation.variable_lower_bound(variable_index);
         const double upper_bound = this->first_reformulation.variable_upper_bound(variable_index);
         if (is_finite(lower_bound)) {
            barrier_terms -= std::log(iterate.primals[variable_index] - lower_bound);
            if (!is_finite(upper_bound)) {
               // damping
               barrier_terms += this->parameters.damping_factor*(iterate.primals[variable_index] - lower_bound);
            }
         }
         if (is_finite(upper_bound)) {
            barrier_terms -= std::log(upper_bound - iterate.primals[variable_index]);
            if (!is_finite(lower_bound)) {
               barrier_terms += this->parameters.damping_factor*(upper_bound - iterate.primals[variable_index]);
            }
         }
      }
      barrier_terms *= this->barrier_parameter;
      assert(!std::isnan(barrier_terms) && "The auxiliary measure is not an number.");
      iterate.progress.auxiliary = barrier_terms;
   }

   double PrimalDualInteriorPointProblem::dual_regularization_factor() const {
      return std::pow(this->barrier_parameter, this->parameters.dual_regularization_exponent);
   }

   // protected member functions

   double PrimalDualInteriorPointProblem::push_variable_to_interior(double variable_value, double lower_bound, double upper_bound) const {
      const double range = upper_bound - lower_bound;
      const double perturbation_lb = std::min(this->parameters.push_variable_to_interior_k1 * std::max(1., std::abs(lower_bound)),
         this->parameters.push_variable_to_interior_k2 * range);
      const double perturbation_ub = std::min(this->parameters.push_variable_to_interior_k1 * std::max(1., std::abs(upper_bound)),
         this->parameters.push_variable_to_interior_k2 * range);
      variable_value = std::max(variable_value, lower_bound + perturbation_lb);
      variable_value = std::min(variable_value, upper_bound - perturbation_ub);
      return variable_value;
   }

   void PrimalDualInteriorPointProblem::compute_bound_dual_direction(const Iterate& current_iterate,
         Direction& direction) const {
      direction.multipliers.lower_bounds.fill(0.);
      direction.multipliers.upper_bounds.fill(0.);
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const double lower_bound = this->first_reformulation.variable_lower_bound(variable_index);
         const double upper_bound = this->first_reformulation.variable_upper_bound(variable_index);
         if (is_finite(lower_bound)) {
            const double distance_to_bound = current_iterate.primals[variable_index] - lower_bound;
            direction.multipliers.lower_bounds[variable_index] = (this->barrier_parameter - direction.primals[variable_index] *
               current_iterate.multipliers.lower_bounds[variable_index]) / distance_to_bound - current_iterate.multipliers.lower_bounds[variable_index];
            assert(is_finite(direction.multipliers.lower_bounds[variable_index]) && "The lower bound dual is infinite");
         }
         if (is_finite(upper_bound)) {
            const double distance_to_bound = current_iterate.primals[variable_index] - upper_bound;
            direction.multipliers.upper_bounds[variable_index] = (this->barrier_parameter - direction.primals[variable_index] *
               current_iterate.multipliers.upper_bounds[variable_index]) / distance_to_bound - current_iterate.multipliers.upper_bounds[variable_index];
            assert(is_finite(direction.multipliers.upper_bounds[variable_index]) && "The upper bound dual is infinite");
         }
      }
   }

   // TODO use a single function for primal and dual fraction-to-boundary rules
   double PrimalDualInteriorPointProblem::primal_fraction_to_boundary(const Vector<double>& current_primals,
         const Vector<double>& primal_direction, double tau) const {
      double step_length = 1.;
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const double lower_bound = this->first_reformulation.variable_lower_bound(variable_index);
         const double upper_bound = this->first_reformulation.variable_upper_bound(variable_index);
         if (is_finite(lower_bound) && primal_direction[variable_index] < 0.) {
            const double distance = -tau * (current_primals[variable_index] - lower_bound) / primal_direction[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
         if (is_finite(upper_bound) && 0. < primal_direction[variable_index]) {
            const double distance = -tau * (current_primals[variable_index] - upper_bound) / primal_direction[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      assert(0. < step_length && step_length <= 1. && "The primal fraction-to-boundary step length is not in (0, 1]");
      return step_length;
   }

   double PrimalDualInteriorPointProblem::dual_fraction_to_boundary(const Multipliers& current_multipliers,
         const Multipliers& direction_multipliers, double tau) const {
      double step_length = 1.;
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const double lower_bound = this->first_reformulation.variable_lower_bound(variable_index);
         const double upper_bound = this->first_reformulation.variable_upper_bound(variable_index);
         if (is_finite(lower_bound) && direction_multipliers.lower_bounds[variable_index] < 0.) {
            const double distance = -tau * current_multipliers.lower_bounds[variable_index] / direction_multipliers.lower_bounds[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
         if (is_finite(upper_bound) && 0. < direction_multipliers.upper_bounds[variable_index]) {
            const double distance = -tau * current_multipliers.upper_bounds[variable_index] / direction_multipliers.upper_bounds[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      assert(0. < step_length && step_length <= 1. && "The dual fraction-to-boundary step length is not in (0, 1]");
      return step_length;
   }

   double PrimalDualInteriorPointProblem::compute_barrier_term_directional_derivative(const Iterate& current_iterate,
         const Vector<double>& primal_direction) const {
      double directional_derivative = 0.;
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const double lower_bound = this->first_reformulation.variable_lower_bound(variable_index);
         const double upper_bound = this->first_reformulation.variable_upper_bound(variable_index);
         if (is_finite(lower_bound)) {
            directional_derivative += -this->barrier_parameter / (current_iterate.primals[variable_index] -
               lower_bound) * primal_direction[variable_index];
            if (!is_finite(upper_bound)) {
               // damping
               directional_derivative += this->parameters.damping_factor * this->barrier_parameter * primal_direction[variable_index];
            }
         }
         if (is_finite(upper_bound)) {
            directional_derivative += -this->barrier_parameter / (current_iterate.primals[variable_index] -
               upper_bound) * primal_direction[variable_index];
            if (!is_finite(lower_bound)) {
               // damping
               directional_derivative -= this->parameters.damping_factor * this->barrier_parameter * primal_direction[variable_index];
            }
         }
      }
      return directional_derivative;
   }

   void PrimalDualInteriorPointProblem::postprocess_iterate(Iterate& iterate) const {
      // rescale the bound multipliers (Eq. 16 in Ipopt paper)
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const double lower_bound = this->first_reformulation.variable_lower_bound(variable_index);
         const double upper_bound = this->first_reformulation.variable_upper_bound(variable_index);
         if (is_finite(lower_bound)) {
            const double coefficient = this->barrier_parameter / (iterate.primals[variable_index] - lower_bound);
            if (is_finite(coefficient)) {
               const double lb = coefficient / this->parameters.k_sigma;
               const double ub = coefficient * this->parameters.k_sigma;
               assert(lb <= ub && "Barrier subproblem: the bounds are in the wrong order in the lower bound multiplier reset");
               if (lb <= ub) {
                  const double current_value = iterate.multipliers.lower_bounds[variable_index];
                  iterate.multipliers.lower_bounds[variable_index] = std::max(std::min(iterate.multipliers.lower_bounds[variable_index], ub), lb);
                  if (iterate.multipliers.lower_bounds[variable_index] != current_value) {
                     DEBUG << "Multiplier for lower bound " << variable_index << " rescaled from " << current_value << " to " <<
                        iterate.multipliers.lower_bounds[variable_index] << '\n';
                  }
               }
               else {
                  WARNING << "Barrier subproblem: the bounds are in the wrong order in the lower bound multiplier reset\n";
               }
            }
         }
         if (is_finite(upper_bound)) {
            const double coefficient = this->barrier_parameter / (iterate.primals[variable_index] - upper_bound);
            if (is_finite(coefficient)) {
               const double lb = coefficient * this->parameters.k_sigma;
               const double ub = coefficient / this->parameters.k_sigma;
               assert(lb <= ub && "Barrier subproblem: the bounds are in the wrong order in the upper bound multiplier reset");
               if (lb <= ub) {
                  const double current_value = iterate.multipliers.upper_bounds[variable_index];
                  iterate.multipliers.upper_bounds[variable_index] = std::max(std::min(iterate.multipliers.upper_bounds[variable_index], ub), lb);
                  if (iterate.multipliers.upper_bounds[variable_index] != current_value) {
                     DEBUG << "Multiplier for upper bound " << variable_index << " rescaled from " << current_value << " to " <<
                        iterate.multipliers.upper_bounds[variable_index] << '\n';
                  }
               }
               else {
                  WARNING << "Barrier subproblem: the bounds are in the wrong order in the upper bound multiplier reset\n";
               }
            }
         }
      }
   }

   double PrimalDualInteriorPointProblem::compute_centrality_error(const Vector<double>& primals,
         const Multipliers& multipliers, double shift) const {
      const Range variables_range = Range(this->first_reformulation.number_variables);
      const VectorExpression shifted_bound_complementarity{variables_range, [&](size_t variable_index) {
         double result = 0.;
         if (0. < multipliers.lower_bounds[variable_index]) { // lower bound
            result = std::max(result, std::abs(multipliers.lower_bounds[variable_index] *
               (primals[variable_index] - this->first_reformulation.variable_lower_bound(variable_index)) - shift));
         }
         if (multipliers.upper_bounds[variable_index] < 0.) { // upper bound
            result = std::max(result, std::abs(multipliers.upper_bounds[variable_index] *
               (primals[variable_index] - this->first_reformulation.variable_upper_bound(variable_index)) - shift));
         }
         return result;
      }};
      return norm_inf(shifted_bound_complementarity); // TODO use a generic norm
   }
} // namespace