// Copyright (c) 2024-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iomanip>
#include "PrimalDualInteriorPointProblem.hpp"
#include "../InteriorPointParameters.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/View.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Parameterization.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   PrimalDualInteriorPointProblem::PrimalDualInteriorPointProblem(const OptimizationProblem& problem,
      const InteriorPointParameters& parameters, const Parameterization& parameterization):
         OptimizationProblem(problem.model, problem.number_variables + problem.get_inequality_constraints().size(),
            problem.number_constraints),
         inner(problem),
         parameterization(parameterization),
         parameters(parameters),
         slacks(this->inner.get_inequality_constraints().size()),
         // all constraints are equality constraints
         equality_constraints(this->number_constraints),
         explicit_variables_lower_bounds(this->number_variables, -INF<double>),
         explicit_variables_upper_bounds(this->number_variables, INF<double>),
         constraints_lower_bounds(this->number_constraints, 0.),
         constraints_upper_bounds(this->number_constraints, 0.),
         variables_lower_bounds(this->number_variables, -INF<double>),
         variables_upper_bounds(this->number_variables, INF<double>) {
      // compute the variables bounds
      view(this->variables_lower_bounds, 0, this->inner.number_variables) = this->inner.get_variables_lower_bounds();
      view(this->variables_upper_bounds, 0, this->inner.number_variables) = this->inner.get_variables_upper_bounds();

      // register the inequality constraint of each slack
      size_t inequality_index = 0;
      for (const size_t constraint_index: this->inner.get_inequality_constraints()) {
         const size_t slack_index = this->inner.number_variables + inequality_index;
         this->slacks.insert(constraint_index, slack_index);
         this->variables_lower_bounds[slack_index] = this->inner.get_constraints_lower_bounds()[constraint_index];
         this->variables_upper_bounds[slack_index] = this->inner.get_constraints_upper_bounds()[constraint_index];
         ++inequality_index;
      }

      // compute the Jacobian sparsity
      const size_t number_jacobian_nonzeros = this->inner.number_jacobian_nonzeros();
      this->jacobian_row_indices.resize(number_jacobian_nonzeros + this->slacks.size());
      this->jacobian_column_indices.resize(number_jacobian_nonzeros + this->slacks.size());
      view(this->jacobian_row_indices, 0, number_jacobian_nonzeros) = this->inner.get_jacobian_row_indices();
      view(this->jacobian_column_indices, 0, number_jacobian_nonzeros) = this->inner.get_jacobian_column_indices();
      size_t nonzero_index = number_jacobian_nonzeros;
      for (const auto [constraint_index, slack_index]: this->slacks) {
         this->jacobian_row_indices[nonzero_index] = static_cast<uno_int>(constraint_index);
         this->jacobian_column_indices[nonzero_index] = static_cast<uno_int>(slack_index);
         ++nonzero_index;
      }
   }

   double PrimalDualInteriorPointProblem::get_objective_multiplier() const {
      return this->inner.get_objective_multiplier();
   }

   bool PrimalDualInteriorPointProblem::has_inequality_constraints() const {
      return false;
   }

   bool PrimalDualInteriorPointProblem::has_bound_constraints() const {
      return false;
   }

   void PrimalDualInteriorPointProblem::generate_initial_iterate(Iterate& initial_iterate, Evaluations& evaluations) const {
      initial_iterate.set_number_variables(this->number_variables);

      // make the initial point strictly feasible wrt the bounds
      for (size_t variable_index: Range(this->inner.number_variables)) {
         initial_iterate.primals[variable_index] = this->push_variable_to_interior(initial_iterate.primals[variable_index],
            this->variables_lower_bounds[variable_index], this->variables_upper_bounds[variable_index]);
      }

      // set the slack variables (if any)
      if (!this->slacks.is_empty()) {
         Vector<double> constraints(this->inner.number_constraints);
         this->inner.evaluate_constraints(initial_iterate, constraints.data(), evaluations);
         // set the slacks to the constraint values
         for (const auto [constraint_index, slack_index]: this->slacks) {
            initial_iterate.primals[slack_index] = this->push_variable_to_interior(constraints[constraint_index],
               this->variables_lower_bounds[slack_index], this->variables_upper_bounds[slack_index]);
         }
         // since the slacks have been set, the function evaluations should also be updated
         evaluations.are_constraints_computed = false;
         evaluations.is_objective_gradient_computed = false;
         evaluations.is_jacobian_computed = false;
      }

      // set the bound multipliers
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            initial_iterate.multipliers.lower_bounds[variable_index] = this->parameters.default_multiplier;
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            initial_iterate.multipliers.upper_bounds[variable_index] = -this->parameters.default_multiplier;
         }
      }
   }

   size_t PrimalDualInteriorPointProblem::number_jacobian_nonzeros() const {
      return this->inner.number_jacobian_nonzeros() + this->slacks.size();
   }

   bool PrimalDualInteriorPointProblem::has_curvature(const HessianModel& hessian_model) const {
      if (hessian_model.has_curvature()) {
         return true;
      }
      else {
         // barrier terms
         for (size_t variable_index: Range(this->number_variables)) {
            if (is_finite(this->variables_lower_bounds[variable_index]) || is_finite(this->variables_upper_bounds[variable_index])) {
               return true;
            }
         }
         return false;
      }
   }

   size_t PrimalDualInteriorPointProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      size_t number_nonzeros = this->inner.number_hessian_nonzeros(hessian_model);
      // barrier contribution
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index]) || is_finite(this->variables_upper_bounds[variable_index])) {
            ++number_nonzeros;
         }
      }
      return number_nonzeros;
   }

   View<const uno_int> PrimalDualInteriorPointProblem::get_jacobian_row_indices() const {
      return this->jacobian_row_indices.view();
   }

   View<const uno_int> PrimalDualInteriorPointProblem::get_jacobian_column_indices() const {
      return this->jacobian_column_indices.view();
   }

   void PrimalDualInteriorPointProblem::compute_hessian_sparsity(const HessianModel& hessian_model, uno_int* row_indices,
         uno_int* column_indices, uno_int solver_indexing) const {
      // original Lagrangian Hessian
      this->inner.compute_hessian_sparsity(hessian_model, row_indices, column_indices, solver_indexing);

      // diagonal barrier terms
      size_t current_index = this->inner.number_hessian_nonzeros(hessian_model);
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index]) || is_finite(this->variables_upper_bounds[variable_index])) {
            row_indices[current_index] = static_cast<uno_int>(variable_index) + solver_indexing;
            column_indices[current_index] = static_cast<uno_int>(variable_index) + solver_indexing;
            ++current_index;
         }
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_constraints(const Iterate& iterate, double* constraints, Evaluations& evaluations) const {
      this->inner.evaluate_constraints(iterate, constraints, evaluations);

      // inequality constraints: add the slacks
      for (const auto [constraint_index, slack_index]: this->slacks) {
         constraints[constraint_index] -= iterate.primals[slack_index];
      }

      // equality constraints: make sure they are homogeneous (c(x) = 0)
      for (const size_t constraint_index: this->inner.get_equality_constraints()) {
         const double fixed_bound = this->inner.get_constraints_lower_bounds()[constraint_index];
         constraints[constraint_index] -= fixed_bound;
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_objective_gradient(const Iterate& iterate, double* objective_gradient,
         Evaluations& evaluations) const {
      this->inner.evaluate_objective_gradient(iterate, objective_gradient, evaluations);

      // barrier terms
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      for (size_t variable_index: Range(this->number_variables)) {
         double barrier_term = 0.;
         if (is_finite(this->variables_lower_bounds[variable_index])) { // lower bounded
            barrier_term += -barrier_parameter/(iterate.primals[variable_index] - this->variables_lower_bounds[variable_index]);
            // damping
            if (is_infinite(this->variables_upper_bounds[variable_index])) {
               barrier_term += this->parameters.damping_factor * barrier_parameter;
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) { // upper bounded
            barrier_term += -barrier_parameter/(iterate.primals[variable_index] - this->variables_upper_bounds[variable_index]);
            // damping
            if (is_infinite(this->variables_lower_bounds[variable_index])) {
               barrier_term -= this->parameters.damping_factor * barrier_parameter;
            }
         }
         objective_gradient[variable_index] += barrier_term;
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_jacobian(const Vector<double>& primals, View<double> jacobian_values,
         Evaluations& evaluations) const {
      this->inner.evaluate_jacobian(primals, jacobian_values, evaluations);

      // add the slack contributions
      size_t nonzero_index = this->inner.number_jacobian_nonzeros();
      for ([[maybe_unused]] const auto _: this->slacks) {
         jacobian_values[nonzero_index] = -1.;
         ++nonzero_index;
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_gradient(const Iterate& iterate, Evaluations& evaluations,
         Vector<double>& lagrangian_gradient) const {
      this->inner.evaluate_lagrangian_gradient(iterate, evaluations, lagrangian_gradient);

      // bound multipliers for slacks
      for (const auto [constraint_index, slack_index]: this->slacks) {
         lagrangian_gradient[slack_index] -= (iterate.multipliers.lower_bounds[slack_index] + iterate.multipliers.upper_bounds[slack_index]);
      }

      // Jacobian block for slacks
      // std::cout << "PrimalDualInteriorPointProblem::evaluate_lagrangian_gradient\n";
      for (const auto [constraint_index, slack_index]: this->slacks) {
         lagrangian_gradient[slack_index] += iterate.multipliers.constraints[constraint_index];
      }

      // barrier terms
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      for (size_t variable_index: Range(this->number_variables)) {
         double barrier_term = 0.;
         if (is_finite(this->variables_lower_bounds[variable_index])) { // lower bounded
            barrier_term += -barrier_parameter/(iterate.primals[variable_index] - this->variables_lower_bounds[variable_index]);
            // damping
            if (is_infinite(this->variables_upper_bounds[variable_index])) {
               barrier_term += this->parameters.damping_factor * barrier_parameter;
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) { // upper bounded
            barrier_term += -barrier_parameter/(iterate.primals[variable_index] - this->variables_upper_bounds[variable_index]);
            // damping
            if (is_infinite(this->variables_lower_bounds[variable_index])) {
               barrier_term -= this->parameters.damping_factor * barrier_parameter;
            }
         }
         // the objective contribution of the Lagrangian gradient may be scaled. Barrier terms go into the constraint contribution
         lagrangian_gradient[variable_index] += barrier_term;
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, View<double> hessian_values) const {
      // original Lagrangian Hessian
      this->inner.evaluate_lagrangian_hessian(statistics, hessian_model, primal_variables, multipliers, hessian_values);

      // barrier terms
      size_t nonzero_index = this->inner.number_hessian_nonzeros(hessian_model);
      for (size_t variable_index: Range(this->number_variables)) {
         const bool finite_lower_bound = is_finite(this->variables_lower_bounds[variable_index]);
         const bool finite_upper_bound = is_finite(this->variables_upper_bounds[variable_index]);
         if (finite_lower_bound || finite_upper_bound) {
            double diagonal_barrier_term = 0.;
            if (finite_lower_bound) {
               const double distance_to_bound = primal_variables[variable_index] - this->variables_lower_bounds[variable_index];
               diagonal_barrier_term += multipliers.lower_bounds[variable_index] / distance_to_bound;
            }
            if (finite_upper_bound) {
               const double distance_to_bound = primal_variables[variable_index] - this->variables_upper_bounds[variable_index];
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

      // add the slack contributions
      for (const auto [constraint_index, slack_variable_index]: this->slacks) {
         result[constraint_index] -= vector[slack_variable_index];
      }
   }

   void PrimalDualInteriorPointProblem::compute_jacobian_transposed_vector_product(const double* vector, double* result,
         const Evaluations& evaluations) const {
      this->inner.compute_jacobian_transposed_vector_product(vector, result, evaluations);

      // add the slack contributions
      for (const auto [constraint_index, slack_variable_index]: this->slacks) {
         result[slack_variable_index] = -vector[constraint_index];
      }
   }

   void PrimalDualInteriorPointProblem::compute_hessian_vector_product(HessianModel& hessian_model, const double* x,
         const double* vector, const Multipliers& multipliers, double* result) const {
      // original Lagrangian Hessian
      this->inner.compute_hessian_vector_product(hessian_model, x, vector, multipliers, result);

      // barrier terms
      for (size_t variable_index: Range(this->number_variables)) {
         const bool finite_lower_bound = is_finite(this->variables_lower_bounds[variable_index]);
         const bool finite_upper_bound = is_finite(this->variables_upper_bounds[variable_index]);
         if (finite_lower_bound || finite_upper_bound) {
            double diagonal_barrier_term = 0.;
            if (finite_lower_bound) { // lower bounded
               const double distance_to_bound = x[variable_index] - this->variables_lower_bounds[variable_index];
               diagonal_barrier_term += multipliers.lower_bounds[variable_index] / distance_to_bound;
            }
            if (finite_upper_bound) { // upper bounded
               const double distance_to_bound = x[variable_index] - this->variables_upper_bounds[variable_index];
               diagonal_barrier_term += multipliers.upper_bounds[variable_index] / distance_to_bound;
            }
            result[variable_index] += diagonal_barrier_term * vector[variable_index];
         }
      }
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_variables_lower_bounds() const {
      return this->explicit_variables_lower_bounds;
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_variables_upper_bounds() const {
      return this->explicit_variables_upper_bounds;
   }

   const Vector<size_t>& PrimalDualInteriorPointProblem::get_fixed_variables() const {
      return this->fixed_variables;
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_constraints_lower_bounds() const {
      return this->constraints_lower_bounds;
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_constraints_upper_bounds() const {
      return this->constraints_upper_bounds;
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_equality_constraints() const {
      return this->equality_constraints;
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_inequality_constraints() const {
      return this->inequality_constraints;
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_dual_regularization_constraints() const {
      return this->inner.get_dual_regularization_constraints();
   }

   Inertia PrimalDualInteriorPointProblem::get_inertia() const {
      return {this->number_variables, this->number_constraints, 0};
   }

   void PrimalDualInteriorPointProblem::assemble_primal_dual_direction(const Iterate& current_iterate, const Vector<double>& solution,
         Direction& direction) const {
      // assemble the primal direction and the constraint dual solution
      OptimizationProblem::assemble_primal_dual_direction(current_iterate, solution, direction);

      // compute the bound duals
      compute_bound_dual_direction(current_iterate, direction);

      // "fraction-to-boundary" rule for primal variables and bound constraints multipliers
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      const double tau = std::max(this->parameters.tau_min, 1. - barrier_parameter);
      direction.primal_dual_step_length = primal_fraction_to_boundary(current_iterate.primals, direction.primals, tau);
      direction.bound_dual_step_length = dual_fraction_to_boundary(current_iterate.multipliers, direction.multipliers, tau);
      DEBUG << "Fraction-to-boundary rules:\n";
      DEBUG << "primal step length = " << direction.primal_dual_step_length << '\n';
      DEBUG << "bound dual step length = " << direction.bound_dual_step_length << "\n\n";
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
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            const double distance_to_bound = current_iterate.primals[variable_index] - this->variables_lower_bounds[variable_index];
            direction.multipliers.lower_bounds[variable_index] = (barrier_parameter - direction.primals[variable_index] *
               current_iterate.multipliers.lower_bounds[variable_index]) / distance_to_bound - current_iterate.multipliers.lower_bounds[variable_index];
            if (is_infinite(direction.multipliers.lower_bounds[variable_index])) {
               throw std::runtime_error("The lower bound dual is infinite");
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            const double distance_to_bound = current_iterate.primals[variable_index] - this->variables_upper_bounds[variable_index];
            direction.multipliers.upper_bounds[variable_index] = (barrier_parameter - direction.primals[variable_index] *
               current_iterate.multipliers.upper_bounds[variable_index]) / distance_to_bound - current_iterate.multipliers.upper_bounds[variable_index];
            if (is_infinite(direction.multipliers.upper_bounds[variable_index])) {
               throw std::runtime_error("The upper bound dual is infinite");
            }
         }
      }
   }

   // TODO use a single function for primal and dual fraction-to-boundary rules
   double PrimalDualInteriorPointProblem::primal_fraction_to_boundary(const Vector<double>& current_primals,
         const Vector<double>& primal_direction, double tau) const {
      // std::cout << "tau = " << tau << '\n';
      double step_length = 1.;
      // original variables
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index]) && primal_direction[variable_index] < 0.) {
            // std::cout << "x_" << variable_index << " has a lower bound " << this->variables_lower_bounds[variable_index] << '\n';
            // std::cout << "d_" << variable_index << " = " << primal_direction[variable_index] << '\n';
            const double distance = -tau * (current_primals[variable_index] - this->variables_lower_bounds[variable_index]) /
               primal_direction[variable_index];
            // std::cout << "For this variable, fraction-to-boundary step length: " << distance << '\n';
            const double new_component = current_primals[variable_index] + distance * primal_direction[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index]) && 0. < primal_direction[variable_index]) {
            const double distance = -tau * (current_primals[variable_index] - this->variables_upper_bounds[variable_index]) /
               primal_direction[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      if (step_length <= 0. || step_length > 1.) {
         throw std::runtime_error("The primal fraction-to-boundary step length is not in (0, 1]");
      }
      // std::cout << "Final fraction-to-boundary step length " << step_length << '\n';
      // std::cout << "New x_0 = " << current_primals[0] + step_length * primal_direction[0] << '\n';
      return step_length;
   }

   double PrimalDualInteriorPointProblem::dual_fraction_to_boundary(const Multipliers& current_multipliers,
         const Multipliers& direction_multipliers, double tau) const {
      double step_length = 1.;
      for (size_t variable_index: Range(this->number_variables)) {
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
      if (step_length <= 0. || step_length > 1.) {
         throw std::runtime_error("The dual fraction-to-boundary step length is not in (0, 1]");
      }
      return step_length;
   }

   double PrimalDualInteriorPointProblem::compute_barrier_term_directional_derivative(const Iterate& current_iterate,
         const Vector<double>& primal_direction) const {
      double directional_derivative = 0.;
      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            directional_derivative += -barrier_parameter / (current_iterate.primals[variable_index] -
               this->variables_lower_bounds[variable_index]) * primal_direction[variable_index];
            if (is_infinite(variables_upper_bounds[variable_index])) {
               // damping
               directional_derivative += this->parameters.damping_factor * barrier_parameter * primal_direction[variable_index];
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            directional_derivative += -barrier_parameter / (current_iterate.primals[variable_index] -
               variables_upper_bounds[variable_index]) * primal_direction[variable_index];
            if (is_infinite(this->variables_lower_bounds[variable_index])) {
               // damping
               directional_derivative -= this->parameters.damping_factor * barrier_parameter * primal_direction[variable_index];
            }
         }
      }
      return directional_derivative;
   }

   void PrimalDualInteriorPointProblem::postprocess_iterate(Iterate& iterate) const {
      const double barrier_parameter = this->parameterization.get("barrier_parameter");

      // if the primals are too close to their bounds, push the bounds away by a small fraction (Section 3.5)
      constexpr double macheps = std::numeric_limits<double>::epsilon();
      const double safe = std::pow(macheps, 0.75);
      size_t adjusted = 0;
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            const double floor = safe * std::max(1.0, std::abs(this->variables_lower_bounds[variable_index]));
            if (iterate.primals[variable_index] - this->variables_lower_bounds[variable_index] < macheps * barrier_parameter) {
               this->variables_lower_bounds[variable_index] = iterate.primals[variable_index] - floor;
               ++adjusted;
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            const double floor = safe * std::max(1.0, std::abs(this->variables_upper_bounds[variable_index]));
            if (this->variables_upper_bounds[variable_index] - iterate.primals[variable_index] < macheps * barrier_parameter) {
               this->variables_upper_bounds[variable_index] = iterate.primals[variable_index] + floor;
               ++adjusted;
            }
         }
      }
      if (adjusted > 0) {
         DEBUG << adjusted << " slack(s) too small, adjusting variable bound\n";
      }

      // rescale the bound multipliers (Eq. 16 in Ipopt paper)
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            const double coefficient = barrier_parameter / (iterate.primals[variable_index] - this->variables_lower_bounds[variable_index]);
            if (is_finite(coefficient)) {
               const double lb = coefficient / this->parameters.k_sigma;
               const double ub = coefficient * this->parameters.k_sigma;
               if (lb > ub) {
                  throw std::runtime_error("Barrier subproblem: the bounds are in the wrong order in the lower bound multiplier reset");
               }
               const double current_value = iterate.multipliers.lower_bounds[variable_index];
               iterate.multipliers.lower_bounds[variable_index] = std::max(std::min(iterate.multipliers.lower_bounds[variable_index], ub), lb);
               if (iterate.multipliers.lower_bounds[variable_index] != current_value) {
                  DEBUG3 << "Multiplier for lower bound " << variable_index << " rescaled from " << current_value << " to " <<
                     iterate.multipliers.lower_bounds[variable_index] << '\n';
               }
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            const double coefficient = barrier_parameter / (iterate.primals[variable_index] - this->variables_upper_bounds[variable_index]);
            if (is_finite(coefficient)) {
               const double lb = coefficient * this->parameters.k_sigma;
               const double ub = coefficient / this->parameters.k_sigma;
               if (lb > ub) {
                  throw std::runtime_error("Barrier subproblem: the bounds are in the wrong order in the upper bound multiplier reset");
               }
               const double current_value = iterate.multipliers.upper_bounds[variable_index];
               iterate.multipliers.upper_bounds[variable_index] = std::max(std::min(iterate.multipliers.upper_bounds[variable_index], ub), lb);
               if (iterate.multipliers.upper_bounds[variable_index] != current_value) {
                  DEBUG3 << "Multiplier for upper bound " << variable_index << " rescaled from " << current_value << " to " <<
                     iterate.multipliers.upper_bounds[variable_index] << '\n';
               }
            }
         }
      }
   }

   static double constraint_violation(const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds,
         double constraint_value, size_t constraint_index) {
      const double lower_bound_violation = std::max(0., lower_bounds[constraint_index] - constraint_value);
      const double upper_bound_violation = std::max(0., constraint_value - upper_bounds[constraint_index]);
      return std::max(lower_bound_violation, upper_bound_violation);
   }

   void PrimalDualInteriorPointProblem::set_infeasibility_measure(Iterate& iterate, Evaluations& evaluations, Norm progress_norm) const {
      Vector<double> constraints(this->number_constraints);
      this->evaluate_constraints(iterate, constraints.data(), evaluations);
      const Range constraints_range = Range(this->number_constraints);
      const VectorExpression v{constraints_range, [&](size_t constraint_index) {
         return constraint_violation(this->constraints_lower_bounds, this->constraints_upper_bounds, constraints[constraint_index],
            constraint_index);
      }};
      iterate.progress.infeasibility = norm(progress_norm, v);
   }

   void PrimalDualInteriorPointProblem::set_objective_measure(Iterate& iterate, Evaluations& evaluations) const {
      this->inner.set_objective_measure(iterate, evaluations);
   }

   void PrimalDualInteriorPointProblem::set_auxiliary_measure(Iterate& iterate) const {
      // start with the auxiliary measure of the initial problem
      this->inner.set_auxiliary_measure(iterate);

      const double barrier_parameter = this->parameterization.get("barrier_parameter");
      if (is_infinite(barrier_parameter)) {
         throw std::runtime_error("Barrier parameter is infinite");
      }

      double barrier_terms = 0.;
      // barrier terms
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            assert(this->variables_lower_bounds[variable_index] < iterate.primals[variable_index]);
            barrier_terms -= std::log(iterate.primals[variable_index] - this->variables_lower_bounds[variable_index]);
            if (is_infinite(this->variables_upper_bounds[variable_index])) {
               // damping
               barrier_terms += this->parameters.damping_factor*(iterate.primals[variable_index] - this->variables_lower_bounds[variable_index]);
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            assert(iterate.primals[variable_index] < this->variables_upper_bounds[variable_index]);
            barrier_terms -= std::log(this->variables_upper_bounds[variable_index] - iterate.primals[variable_index]);
            if (is_infinite(this->variables_lower_bounds[variable_index])) {
               // damping
               barrier_terms += this->parameters.damping_factor*(this->variables_upper_bounds[variable_index] - iterate.primals[variable_index]);
            }
         }
      }
      barrier_terms *= barrier_parameter;
      if (std::isnan(barrier_terms)) {
         throw std::runtime_error("The auxiliary measure is not an number.");
      }
      iterate.progress.auxiliary += barrier_terms;
   }

   // predicted reductions

   double PrimalDualInteriorPointProblem::compute_predicted_infeasibility_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length, Norm norm, Evaluations& current_evaluations) const {
      return this->inner.compute_predicted_infeasibility_reduction(current_iterate, primal_direction, step_length, norm,
         current_evaluations);
   }

   std::function<double(double)> PrimalDualInteriorPointProblem::compute_predicted_objective_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length, Evaluations& current_evaluations,
         double hessian_quadratic_form) const {
      return this->inner.compute_predicted_objective_reduction(current_iterate, primal_direction, step_length,
         current_evaluations, hessian_quadratic_form);
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

   double PrimalDualInteriorPointProblem::compute_centrality_error(const Vector<double>& primals, const Multipliers& multipliers,
         double shift) const {
      const Range variables_range = Range(this->number_variables);
      const VectorExpression shifted_bound_complementarity{variables_range, [&](size_t variable_index) {
         double result = 0.;
         if (0. < multipliers.lower_bounds[variable_index]) { // lower bound
            result = std::max(result, std::abs(multipliers.lower_bounds[variable_index] *
               (primals[variable_index] - this->variables_lower_bounds[variable_index]) - shift));
         }
         if (multipliers.upper_bounds[variable_index] < 0.) { // upper bound
            result = std::max(result, std::abs(multipliers.upper_bounds[variable_index] *
               (primals[variable_index] - this->variables_upper_bounds[variable_index]) - shift));
         }
         return result;
      }};
      return norm_inf(shifted_bound_complementarity); // TODO use a generic norm
   }
} // namespace