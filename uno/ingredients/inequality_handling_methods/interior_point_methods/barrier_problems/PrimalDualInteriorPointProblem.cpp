// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalDualInteriorPointProblem.hpp"
#include "../InteriorPointParameters.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/VectorView.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Parameterization.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   PrimalDualInteriorPointProblem::PrimalDualInteriorPointProblem(const OptimizationProblem& problem,
      const InteriorPointParameters& parameters, const Parameterization& parameterization):
         OptimizationProblem(problem.model, problem.number_variables, problem.number_constraints),
         inner(problem),
         parameterization(parameterization),
         parameters(parameters),
         equality_constraints(problem.number_constraints),
         barrier_variables_lower_bounds(this->number_variables, -INF<double>),
         barrier_variables_upper_bounds(this->number_variables, INF<double>) {
   }

   std::unique_ptr<OptimizationProblem> PrimalDualInteriorPointProblem::clone() const {
      return std::make_unique<PrimalDualInteriorPointProblem>(*this);
   }

   double PrimalDualInteriorPointProblem::get_objective_multiplier() const {
      return this->inner.get_objective_multiplier();
   }

   void PrimalDualInteriorPointProblem::generate_initial_iterate(Iterate& initial_iterate, Evaluations& evaluations) const {
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();

      // make the initial point strictly feasible wrt the bounds
      for (size_t variable_index: Range(this->number_variables)) {
         initial_iterate.primals[variable_index] = this->push_variable_to_interior(initial_iterate.primals[variable_index],
            variables_lower_bounds[variable_index], variables_upper_bounds[variable_index]);
      }

      // set the slack variables (if any)
      if (!this->model.get_slacks().is_empty()) {
         evaluations.evaluate_constraints(this->model, initial_iterate.primals);
         // set the slacks to the constraint values
         for (const auto [constraint_index, slack_index]: this->model.get_slacks()) {
            initial_iterate.primals[slack_index] = this->push_variable_to_interior(evaluations.constraints[constraint_index],
               variables_lower_bounds[slack_index], variables_upper_bounds[slack_index]);
         }
         // since the slacks have been set, the function evaluations should also be updated
         evaluations.are_constraints_computed = false;
         evaluations.is_objective_gradient_computed = false;
         evaluations.is_jacobian_computed = false;
      }

      // set the bound multipliers
      for (size_t variable_index: Range(this->inner.number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index])) {
            initial_iterate.multipliers.lower_bounds[variable_index] = this->parameters.default_multiplier;
         }
         if (is_finite(variables_upper_bounds[variable_index])) {
            initial_iterate.multipliers.upper_bounds[variable_index] = -this->parameters.default_multiplier;
         }
      }

      // TODO compute least-square multipliers
      if (0 < this->number_constraints) {
      }
   }

   size_t PrimalDualInteriorPointProblem::number_jacobian_nonzeros() const {
      return this->inner.number_jacobian_nonzeros();
   }

   bool PrimalDualInteriorPointProblem::has_curvature(const HessianModel& hessian_model) const {
      if (hessian_model.has_curvature()) {
         return true;
      }
      else {
         // barrier terms
         const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
         const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
         for (size_t variable_index: Range(this->inner.number_variables)) {
            if (is_finite(variables_lower_bounds[variable_index]) || is_finite(variables_upper_bounds[variable_index])) {
               return true;
            }
         }
         return false;
      }
   }

   size_t PrimalDualInteriorPointProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      size_t number_nonzeros = this->inner.number_hessian_nonzeros(hessian_model);
      // barrier contribution: original variables
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index]) || is_finite(variables_upper_bounds[variable_index])) {
            ++number_nonzeros;
         }
      }
      return number_nonzeros;
   }

   void PrimalDualInteriorPointProblem::compute_jacobian_sparsity(uno_int* row_indices, uno_int* column_indices,
         uno_int row_offset, uno_int column_offset, uno_int solver_indexing, MatrixOrder matrix_order) const {
      this->inner.compute_jacobian_sparsity(row_indices, column_indices, row_offset, column_offset,
         solver_indexing, matrix_order);
   }

   void PrimalDualInteriorPointProblem::compute_hessian_sparsity(const HessianModel& hessian_model, uno_int* row_indices,
         uno_int* column_indices, uno_int solver_indexing) const {
      // original Lagrangian Hessian
      this->inner.compute_hessian_sparsity(hessian_model, row_indices, column_indices, solver_indexing);

      // diagonal barrier terms
      size_t current_index = this->inner.number_hessian_nonzeros(hessian_model);
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index]) || is_finite(variables_upper_bounds[variable_index])) {
            row_indices[current_index] = static_cast<int>(variable_index) + solver_indexing;
            column_indices[current_index] = static_cast<int>(variable_index) + solver_indexing;
            ++current_index;
         }
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_constraints(const Iterate& iterate, double* constraints, Evaluations& evaluations) const {
      this->inner.evaluate_constraints(iterate, constraints, evaluations);
   }

   void PrimalDualInteriorPointProblem::evaluate_objective_gradient(const Iterate& iterate, double* objective_gradient,
         Evaluations& evaluations) const {
      this->inner.evaluate_objective_gradient(iterate, objective_gradient, evaluations);

      // barrier terms
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         double barrier_term = 0.;
         if (is_finite(variables_lower_bounds[variable_index])) { // lower bounded
            barrier_term += -barrier_parameter/(iterate.primals[variable_index] - variables_lower_bounds[variable_index]);
            // damping
            if (is_infinite(variables_upper_bounds[variable_index])) {
               barrier_term += this->parameters.damping_factor * barrier_parameter;
            }
         }
         if (is_finite(variables_upper_bounds[variable_index])) { // upper bounded
            barrier_term += -barrier_parameter/(iterate.primals[variable_index] - variables_upper_bounds[variable_index]);
            // damping
            if (is_infinite(variables_lower_bounds[variable_index])) {
               barrier_term -= this->parameters.damping_factor * barrier_parameter;
            }
         }
         objective_gradient[variable_index] += barrier_term;
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_jacobian(const Vector<double>& primals, double* jacobian_values, Evaluations& evaluations) const {
      this->inner.evaluate_jacobian(primals, jacobian_values, evaluations);
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_gradient(const Iterate& iterate, Evaluations& evaluations,
         Vector<double>& lagrangian_gradient) const {
      this->inner.evaluate_lagrangian_gradient(iterate, evaluations, lagrangian_gradient);

      // barrier terms
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         double barrier_term = 0.;
         if (is_finite(variables_lower_bounds[variable_index])) { // lower bounded
            barrier_term += -barrier_parameter/(iterate.primals[variable_index] - variables_lower_bounds[variable_index]);
            // damping
            if (is_infinite(variables_upper_bounds[variable_index])) {
               barrier_term += this->parameters.damping_factor * barrier_parameter;
            }
         }
         if (is_finite(variables_upper_bounds[variable_index])) { // upper bounded
            barrier_term += -barrier_parameter/(iterate.primals[variable_index] - variables_upper_bounds[variable_index]);
            // damping
            if (is_infinite(variables_lower_bounds[variable_index])) {
               barrier_term -= this->parameters.damping_factor * barrier_parameter;
            }
         }
         // the objective contribution of the Lagrangian gradient may be scaled. Barrier terms go into the constraint contribution
         lagrangian_gradient[variable_index] += barrier_term;
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, double* hessian_values) const {
      // original Lagrangian Hessian
      this->inner.evaluate_lagrangian_hessian(statistics, hessian_model, primal_variables, multipliers, hessian_values);

      // barrier terms
      size_t nonzero_index = this->inner.number_hessian_nonzeros(hessian_model);
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         const bool finite_lower_bound = is_finite(variables_lower_bounds[variable_index]);
         const bool finite_upper_bound = is_finite(variables_upper_bounds[variable_index]);
         if (finite_lower_bound || finite_upper_bound) {
            double diagonal_barrier_term = 0.;
            if (finite_lower_bound) {
               const double distance_to_bound = primal_variables[variable_index] - variables_lower_bounds[variable_index];
               diagonal_barrier_term += multipliers.lower_bounds[variable_index] / distance_to_bound;
            }
            if (finite_upper_bound) {
               const double distance_to_bound = primal_variables[variable_index] - variables_upper_bounds[variable_index];
               diagonal_barrier_term += multipliers.upper_bounds[variable_index] / distance_to_bound;
            }
            hessian_values[nonzero_index] = diagonal_barrier_term;
            ++nonzero_index;
         }
      }
   }

   void PrimalDualInteriorPointProblem::compute_jacobian_vector_product(const double* vector, double* result,
         const Evaluations& evaluations) const {
      this->inner.compute_jacobian_vector_product(vector, result, evaluations);
   }

   void PrimalDualInteriorPointProblem::compute_jacobian_transposed_vector_product(const double* vector, double* result,
         const Evaluations& evaluations) const {
      this->inner.compute_jacobian_transposed_vector_product(vector, result, evaluations);
   }

   void PrimalDualInteriorPointProblem::compute_hessian_vector_product(HessianModel& hessian_model, const double* x,
         const double* vector, const Multipliers& multipliers, double* result) const {
      // original Lagrangian Hessian
      this->inner.compute_hessian_vector_product(hessian_model, x, vector, multipliers, result);

      // barrier terms
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         const bool finite_lower_bound = is_finite(variables_lower_bounds[variable_index]);
         const bool finite_upper_bound = is_finite(variables_upper_bounds[variable_index]);
         if (finite_lower_bound || finite_upper_bound) {
            double diagonal_barrier_term = 0.;
            if (finite_lower_bound) { // lower bounded
               const double distance_to_bound = vector[variable_index] - variables_lower_bounds[variable_index];
               diagonal_barrier_term += multipliers.lower_bounds[variable_index] / distance_to_bound;
            }
            if (finite_upper_bound) { // upper bounded
               const double distance_to_bound = vector[variable_index] - variables_upper_bounds[variable_index];
               diagonal_barrier_term += multipliers.upper_bounds[variable_index] / distance_to_bound;
            }
            result[variable_index] += diagonal_barrier_term * vector[variable_index];
         }
      }
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_variables_lower_bounds() const {
      return this->barrier_variables_lower_bounds;
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_variables_upper_bounds() const {
      return this->barrier_variables_lower_bounds;
   }

   const Vector<size_t>& PrimalDualInteriorPointProblem::get_fixed_variables() const {
      return this->fixed_variables;
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_constraints_lower_bounds() const {
      return this->inner.get_constraints_lower_bounds();
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_constraints_upper_bounds() const {
      return this->inner.get_constraints_upper_bounds();
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_equality_constraints() const {
      return this->equality_constraints;
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_inequality_constraints() const {
      return this->inequality_constraints;
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_dual_regularization_constraints() const {
      if (this->inner.get_dual_regularization_constraints().empty()) {
         // this is an indication that the constraints (if there is any) were already regularized in a previous
         // reformulation (e.g. l1 relaxation). In that case, we stick to an empty set
         return this->inner.get_dual_regularization_constraints();
      }
      // otherwise, we pick the set of equality constraints, since the inequality constraints have slacks
      return this->inner.get_equality_constraints();
   }

   void PrimalDualInteriorPointProblem::assemble_primal_dual_direction(const Iterate& current_iterate, const Vector<double>& solution,
         Direction& direction) const {
      // assemble the primal direction and the constraint dual solution
      OptimizationProblem::assemble_primal_dual_direction(current_iterate, solution, direction);

      // compute the bound duals
      this->compute_bound_dual_direction(current_iterate, direction);

      // "fraction-to-boundary" rule for primal variables and constraints multipliers
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      const double tau = std::max(this->parameters.tau_min, 1. - barrier_parameter);
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

   double PrimalDualInteriorPointProblem::dual_regularization_factor() const {
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      return std::pow(barrier_parameter, this->parameters.dual_regularization_exponent);
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
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index])) {
            const double distance_to_bound = current_iterate.primals[variable_index] - variables_lower_bounds[variable_index];
            direction.multipliers.lower_bounds[variable_index] = (barrier_parameter - direction.primals[variable_index] *
               current_iterate.multipliers.lower_bounds[variable_index]) / distance_to_bound - current_iterate.multipliers.lower_bounds[variable_index];
            assert(is_finite(direction.multipliers.lower_bounds[variable_index]) && "The lower bound dual is infinite");
         }
         if (is_finite(variables_upper_bounds[variable_index])) {
            const double distance_to_bound = current_iterate.primals[variable_index] - variables_upper_bounds[variable_index];
            direction.multipliers.upper_bounds[variable_index] = (barrier_parameter - direction.primals[variable_index] *
               current_iterate.multipliers.upper_bounds[variable_index]) / distance_to_bound - current_iterate.multipliers.upper_bounds[variable_index];
            assert(is_finite(direction.multipliers.upper_bounds[variable_index]) && "The upper bound dual is infinite");
         }
      }
   }

   // TODO use a single function for primal and dual fraction-to-boundary rules
   double PrimalDualInteriorPointProblem::primal_fraction_to_boundary(const Vector<double>& current_primals,
         const Vector<double>& primal_direction, double tau) const {
      double step_length = 1.;
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index]) && primal_direction[variable_index] < 0.) {
            const double distance = -tau * (current_primals[variable_index] - variables_lower_bounds[variable_index]) /
               primal_direction[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
         if (is_finite(variables_upper_bounds[variable_index]) && 0. < primal_direction[variable_index]) {
            const double distance = -tau * (current_primals[variable_index] - variables_upper_bounds[variable_index]) /
               primal_direction[variable_index];
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
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index]) && direction_multipliers.lower_bounds[variable_index] < 0.) {
            const double distance = -tau * current_multipliers.lower_bounds[variable_index] / direction_multipliers.lower_bounds[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
         if (is_finite(variables_upper_bounds[variable_index]) && 0. < direction_multipliers.upper_bounds[variable_index]) {
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
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index])) {
            directional_derivative += -barrier_parameter / (current_iterate.primals[variable_index] -
               variables_lower_bounds[variable_index]) * primal_direction[variable_index];
            if (is_infinite(variables_upper_bounds[variable_index])) {
               // damping
               directional_derivative += this->parameters.damping_factor * barrier_parameter * primal_direction[variable_index];
            }
         }
         if (is_finite(variables_upper_bounds[variable_index])) {
            directional_derivative += -barrier_parameter / (current_iterate.primals[variable_index] -
               variables_upper_bounds[variable_index]) * primal_direction[variable_index];
            if (is_infinite(variables_lower_bounds[variable_index])) {
               // damping
               directional_derivative -= this->parameters.damping_factor * barrier_parameter * primal_direction[variable_index];
            }
         }
      }
      return directional_derivative;
   }

   void PrimalDualInteriorPointProblem::postprocess_iterate(Iterate& iterate) const {
      // rescale the bound multipliers (Eq. 16 in Ipopt paper)
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index])) {
            const double coefficient = barrier_parameter / (iterate.primals[variable_index] - variables_lower_bounds[variable_index]);
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
         if (is_finite(variables_upper_bounds[variable_index])) {
            const double coefficient = barrier_parameter / (iterate.primals[variable_index] - variables_upper_bounds[variable_index]);
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

   void PrimalDualInteriorPointProblem::set_auxiliary_measure(Iterate& iterate) const {
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      assert(is_finite(barrier_parameter));

      // start with the auxiliary measure of the initial problem
      this->inner.set_auxiliary_measure(iterate);

      // add the contribution of the barrier terms
      double barrier_terms = 0.;
      const auto& variables_lower_bounds = this->inner.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->inner.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->inner.number_variables)) {
         if (is_finite(variables_lower_bounds[variable_index])) {
            barrier_terms -= std::log(iterate.primals[variable_index] - variables_lower_bounds[variable_index]);
            if (is_infinite(variables_upper_bounds[variable_index])) {
               // damping
               barrier_terms += this->parameters.damping_factor*(iterate.primals[variable_index] - variables_lower_bounds[variable_index]);
            }
         }
         if (is_finite(variables_upper_bounds[variable_index])) {
            barrier_terms -= std::log(variables_upper_bounds[variable_index] - iterate.primals[variable_index]);
            if (is_infinite(variables_lower_bounds[variable_index])) {
               barrier_terms += this->parameters.damping_factor*(variables_upper_bounds[variable_index] - iterate.primals[variable_index]);
            }
         }
      }
      barrier_terms *= barrier_parameter;
      assert(!std::isnan(barrier_terms) && "The auxiliary measure is not an number.");
      iterate.progress.auxiliary += barrier_terms;
   }

   double PrimalDualInteriorPointProblem::compute_predicted_auxiliary_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const {
      // start with the auxiliary measure of the initial problem
      double predicted_auxiliary_reduction = this->inner.compute_predicted_auxiliary_reduction(current_iterate,
         primal_direction, step_length);

      // add the contribution of the barrier terms
      const double directional_derivative = this->compute_barrier_term_directional_derivative(current_iterate, primal_direction);
      predicted_auxiliary_reduction += step_length * (-directional_derivative);
      // }, "α*(μ*X^{-1} e^T d)"};
      return predicted_auxiliary_reduction;
   }
} // namespace