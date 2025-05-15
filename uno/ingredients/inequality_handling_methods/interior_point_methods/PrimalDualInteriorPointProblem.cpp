// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalDualInteriorPointProblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   PrimalDualInteriorPointProblem::PrimalDualInteriorPointProblem(const OptimizationProblem& first_reformulation, double barrier_parameter):
         OptimizationProblem(first_reformulation.model,
            first_reformulation.number_variables + first_reformulation.get_inequality_constraints().size(),
            first_reformulation.number_constraints),
         first_reformulation(first_reformulation),
         number_slack_variables(first_reformulation.get_inequality_constraints().size()),
         equality_constraints(first_reformulation.number_constraints),
         barrier_parameter(barrier_parameter) {
      if (!first_reformulation.get_fixed_variables().empty()) {
         throw std::runtime_error("The problem has fixed variables. Move them to the set of general constraints.");
      }
   }

   double PrimalDualInteriorPointProblem::get_objective_multiplier() const {
      return this->first_reformulation.get_objective_multiplier();
   }

   void PrimalDualInteriorPointProblem::evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const {
      this->first_reformulation.evaluate_objective_gradient(iterate, objective_gradient);

      // barrier terms of original variables
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const bool lower_bounded = is_finite(this->original_variable_lower_bound(variable_index));
         const bool upper_bounded = is_finite(this->original_variable_upper_bound(variable_index));
         if (lower_bounded || upper_bounded) {
            double barrier_term = 0.;
            if (lower_bounded) {
               barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - this->original_variable_lower_bound(variable_index));
               // damping
               if (!upper_bounded) {
                  barrier_term += this->damping_factor * this->barrier_parameter;
               }
            }
            if (upper_bounded) {
               barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - this->original_variable_upper_bound(variable_index));
               // damping
               if (!lower_bounded) {
                  barrier_term -= this->damping_factor * this->barrier_parameter;
               }
            }
            objective_gradient.insert(variable_index, barrier_term);
         }
      }

      // barrier terms of slack variables
      size_t slack_index = this->first_reformulation.number_variables;
      for (size_t inequality_constraint_index: this->first_reformulation.get_inequality_constraints()) {
         const bool lower_bounded = is_finite(this->slack_lower_bound(inequality_constraint_index));
         const bool upper_bounded = is_finite(this->slack_upper_bound(inequality_constraint_index));
         double barrier_term = 0.;
         if (lower_bounded) {
            barrier_term += -this->barrier_parameter/(iterate.primals[slack_index] - this->slack_lower_bound(inequality_constraint_index));
            // damping
            if (!upper_bounded) {
               barrier_term += this->damping_factor * this->barrier_parameter;
            }
         }
         if (upper_bounded) {
            barrier_term += -this->barrier_parameter/(iterate.primals[slack_index] - this->slack_upper_bound(inequality_constraint_index));
            // damping
            if (!lower_bounded) {
               barrier_term -= this->damping_factor * this->barrier_parameter;
            }
         }
         objective_gradient.insert(slack_index, barrier_term);
         slack_index++;
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
      this->first_reformulation.evaluate_constraints(iterate, constraints);
      // inequality constraints l <= c(x) <= u: add the value of the slack. This results in c(x) - s = 0
      size_t slack_index = this->first_reformulation.number_variables;
      for (size_t inequality_constraint_index: this->first_reformulation.get_inequality_constraints()) {
         constraints[inequality_constraint_index] -= iterate.primals[slack_index];
         slack_index++;
      }
      // equality constraints c(x) = l: make sure they are homogeneous (c(x) - l = 0)
      for (const size_t constraint_index: this->first_reformulation.get_equality_constraints()) {
         constraints[constraint_index] -= this->first_reformulation.constraint_lower_bound(constraint_index);
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
      this->first_reformulation.evaluate_constraint_jacobian(iterate, constraint_jacobian);
      // add slack contribution
      size_t slack_index = this->first_reformulation.number_variables;
      for (size_t inequality_constraint_index: this->first_reformulation.get_inequality_constraints()) {
         constraint_jacobian[inequality_constraint_index].insert(slack_index, -1.);
         slack_index++;
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, SymmetricMatrix<size_t, double>& hessian) const {
      // original Lagrangian Hessian
      this->first_reformulation.evaluate_lagrangian_hessian(statistics, hessian_model, primal_variables, multipliers, hessian);
      hessian.set_dimension(this->number_variables);

      // barrier terms of original variables
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const bool lower_bounded = is_finite(this->original_variable_lower_bound(variable_index));
         const bool upper_bounded = is_finite(this->original_variable_upper_bound(variable_index));
         if (lower_bounded || upper_bounded) {
            double diagonal_barrier_term = 0.;
            if (lower_bounded) {
               const double distance_to_bound = primal_variables[variable_index] - this->original_variable_lower_bound(variable_index);
               diagonal_barrier_term += multipliers.lower_bounds[variable_index] / distance_to_bound;
            }
            if (upper_bounded) {
               const double distance_to_bound = primal_variables[variable_index] - this->original_variable_upper_bound(variable_index);
               diagonal_barrier_term += multipliers.upper_bounds[variable_index] / distance_to_bound;
            }
            hessian.insert(diagonal_barrier_term, variable_index, variable_index);
         }
      }

      // barrier terms of slack variables
      size_t slack_index = this->first_reformulation.number_variables;
      for (size_t inequality_constraint_index: this->first_reformulation.get_inequality_constraints()) {
         double diagonal_barrier_term = 0.;
         if (is_finite(this->slack_lower_bound(inequality_constraint_index))) {
            const double distance_to_bound = primal_variables[slack_index] - this->slack_lower_bound(inequality_constraint_index);
            diagonal_barrier_term += multipliers.lower_bounds[slack_index] / distance_to_bound;
         }
         if (is_finite(this->slack_upper_bound(inequality_constraint_index))) {
            const double distance_to_bound = primal_variables[slack_index] - this->slack_upper_bound(inequality_constraint_index);
            diagonal_barrier_term += multipliers.upper_bounds[slack_index] / distance_to_bound;
         }
         hessian.insert(diagonal_barrier_term, slack_index, slack_index);
         slack_index++;
      }
   }

   void PrimalDualInteriorPointProblem::compute_hessian_vector_product(HessianModel& hessian_model, const Vector<double>& vector,
         const Multipliers& multipliers, Vector<double>& result) const {
      // original Lagrangian Hessian
      this->first_reformulation.compute_hessian_vector_product(hessian_model, vector, multipliers, result);

      // barrier terms of original variables
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const bool lower_bounded = is_finite(this->original_variable_lower_bound(variable_index));
         const bool upper_bounded = is_finite(this->original_variable_upper_bound(variable_index));
         if (lower_bounded || upper_bounded) {
            double diagonal_barrier_term = 0.;
            if (lower_bounded) {
               const double distance_to_bound = vector[variable_index] - this->original_variable_lower_bound(variable_index);
               diagonal_barrier_term += multipliers.lower_bounds[variable_index] / distance_to_bound;
            }
            if (upper_bounded) {
               const double distance_to_bound = vector[variable_index] - this->original_variable_upper_bound(variable_index);
               diagonal_barrier_term += multipliers.upper_bounds[variable_index] / distance_to_bound;
            }
            result[variable_index] += diagonal_barrier_term * vector[variable_index];
         }
      }

      // barrier terms of slack variables
      size_t slack_index = this->first_reformulation.number_variables;
      for (size_t inequality_constraint_index: this->first_reformulation.get_inequality_constraints()) {
         double diagonal_barrier_term = 0.;
         if (is_finite(this->slack_lower_bound(inequality_constraint_index))) { // lower bounded
            const double distance_to_bound = vector[slack_index] - this->slack_lower_bound(inequality_constraint_index);
            diagonal_barrier_term += multipliers.lower_bounds[slack_index] / distance_to_bound;
         }
         if (is_finite(this->slack_upper_bound(inequality_constraint_index))) { // upper bounded
            const double distance_to_bound = vector[slack_index] - this->slack_upper_bound(inequality_constraint_index);
            diagonal_barrier_term += multipliers.upper_bounds[slack_index] / distance_to_bound;
         }
         result[slack_index] += diagonal_barrier_term * vector[slack_index];
         slack_index++;
      }
   }

   double PrimalDualInteriorPointProblem::primal_fraction_to_boundary(const Vector<double>& current_primals,
         const Vector<double>& primal_direction, double tau) const {
      double step_length = 1.;
      // original variables
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         if (is_finite(this->original_variable_lower_bound(variable_index)) && primal_direction[variable_index] < 0.) {
            double distance = -tau * (current_primals[variable_index] - this->original_variable_lower_bound(variable_index)) / primal_direction[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
         if (is_finite(this->original_variable_upper_bound(variable_index)) && 0. < primal_direction[variable_index]) {
            double distance = -tau * (current_primals[variable_index] - this->original_variable_upper_bound(variable_index)) / primal_direction[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      // slack variables
      size_t slack_index = this->first_reformulation.number_variables;
      for (size_t inequality_constraint_index: this->first_reformulation.get_inequality_constraints()) {
         if (is_finite(this->slack_lower_bound(inequality_constraint_index)) && primal_direction[slack_index] < 0.) {
            double distance = -tau * (current_primals[slack_index] - this->slack_lower_bound(inequality_constraint_index)) / primal_direction[slack_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
         if (is_finite(this->slack_upper_bound(inequality_constraint_index)) && 0. < primal_direction[slack_index]) {
            double distance = -tau * (current_primals[slack_index] - this->slack_upper_bound(inequality_constraint_index)) / primal_direction[slack_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
         slack_index++;
      }
      assert(0. < step_length && step_length <= 1. && "The primal fraction-to-boundary step length is not in (0, 1]");
      return step_length;
   }

   double PrimalDualInteriorPointProblem::dual_fraction_to_boundary(const Multipliers& current_multipliers,
         Multipliers& direction_multipliers, double tau) const {
      double step_length = 1.;
      // original variables
      for (const size_t variable_index: this->first_reformulation.get_lower_bounded_variables()) {
         if (direction_multipliers.lower_bounds[variable_index] < 0.) {
            double distance = -tau * current_multipliers.lower_bounds[variable_index] / direction_multipliers.lower_bounds[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      for (const size_t variable_index: this->first_reformulation.get_upper_bounded_variables()) {
         if (0. < direction_multipliers.upper_bounds[variable_index]) {
            double distance = -tau * current_multipliers.upper_bounds[variable_index] / direction_multipliers.upper_bounds[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      // slack variables
      size_t slack_index = this->first_reformulation.number_variables;
      for (size_t inequality_constraint_index: this->first_reformulation.get_inequality_constraints()) {
         if (is_finite(this->slack_lower_bound(inequality_constraint_index)) && direction_multipliers.lower_bounds[slack_index] < 0.) {
            double distance = -tau * current_multipliers.lower_bounds[slack_index] / direction_multipliers.lower_bounds[slack_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
         if (is_finite(this->slack_upper_bound(inequality_constraint_index)) && 0. < direction_multipliers.upper_bounds[slack_index]) {
            double distance = -tau * current_multipliers.upper_bounds[slack_index] / direction_multipliers.upper_bounds[slack_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
         slack_index++;
      }
      assert(0. < step_length && step_length <= 1. && "The dual fraction-to-boundary step length is not in (0, 1]");
      return step_length;
   }

   void PrimalDualInteriorPointProblem::compute_bound_dual_direction(const Vector<double>& current_primals,
         const Multipliers& current_multipliers, const Vector<double>& primal_direction, Multipliers& direction_multipliers) const {
      direction_multipliers.lower_bounds.fill(0.);
      direction_multipliers.upper_bounds.fill(0.);
      // original variables
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         if (is_finite(this->original_variable_lower_bound(variable_index))) {
            const double distance_to_bound = current_primals[variable_index] - this->original_variable_lower_bound(variable_index);
            direction_multipliers.lower_bounds[variable_index] = (this->barrier_parameter - primal_direction[variable_index] * current_multipliers.lower_bounds[variable_index]) /
               distance_to_bound - current_multipliers.lower_bounds[variable_index];
            assert(is_finite(direction_multipliers.lower_bounds[variable_index]) && "The lower bound dual is infinite");
         }
         if (is_finite(this->original_variable_upper_bound(variable_index))) {
            const double distance_to_bound = current_primals[variable_index] - this->original_variable_upper_bound(variable_index);
            direction_multipliers.upper_bounds[variable_index] = (this->barrier_parameter - primal_direction[variable_index] * current_multipliers.upper_bounds[variable_index]) /
               distance_to_bound - current_multipliers.upper_bounds[variable_index];
            assert(is_finite(direction_multipliers.upper_bounds[variable_index]) && "The upper bound dual is infinite");
         }
      }
      // slack variables
      size_t slack_index = this->first_reformulation.number_variables;
      for (size_t inequality_constraint_index: this->first_reformulation.get_inequality_constraints()) {
         if (is_finite(this->slack_lower_bound(inequality_constraint_index))) {
            const double distance_to_bound = current_primals[slack_index] - this->slack_lower_bound(inequality_constraint_index);
            direction_multipliers.lower_bounds[slack_index] = (this->barrier_parameter - primal_direction[slack_index] * current_multipliers.lower_bounds[slack_index]) /
               distance_to_bound - current_multipliers.lower_bounds[slack_index];
            assert(is_finite(direction_multipliers.lower_bounds[slack_index]) && "The lower bound dual is infinite");
         }
         if (is_finite(this->slack_upper_bound(inequality_constraint_index))) {
            const double distance_to_bound = current_primals[slack_index] - this->slack_upper_bound(inequality_constraint_index);
            direction_multipliers.upper_bounds[slack_index] = (this->barrier_parameter - primal_direction[slack_index] * current_multipliers.upper_bounds[slack_index]) /
               distance_to_bound - current_multipliers.upper_bounds[slack_index];
            assert(is_finite(direction_multipliers.upper_bounds[slack_index]) && "The upper bound dual is infinite");
         }
         slack_index++;
      }
   }

   double PrimalDualInteriorPointProblem::compute_auxiliary_measure(Iterate& iterate) const {
      // auxiliary measure: barrier terms
      double barrier_terms = 0.;
      // original variables
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const bool lower_bounded = is_finite(this->original_variable_lower_bound(variable_index));
         const bool upper_bounded = is_finite(this->original_variable_upper_bound(variable_index));
         if (lower_bounded) {
            barrier_terms -= std::log(iterate.primals[variable_index] - this->original_variable_lower_bound(variable_index));
            // damping
            if (!upper_bounded) {
               barrier_terms += this->damping_factor*(iterate.primals[variable_index] - this->original_variable_lower_bound(variable_index));
            }
         }
         if (upper_bounded) {
            barrier_terms -= std::log(this->original_variable_upper_bound(variable_index) - iterate.primals[variable_index]);
            // damping
            if (!lower_bounded) {
               barrier_terms += this->damping_factor*(this->original_variable_upper_bound(variable_index) - iterate.primals[variable_index]);
            }
         }
      }
      // slack variables
      size_t slack_index = this->first_reformulation.number_variables;
      for (size_t inequality_constraint_index: this->first_reformulation.get_inequality_constraints()) {
         const bool lower_bounded = is_finite(this->slack_lower_bound(inequality_constraint_index));
         const bool upper_bounded = is_finite(this->slack_upper_bound(inequality_constraint_index));
         if (lower_bounded) {
            barrier_terms -= std::log(iterate.primals[slack_index] - this->slack_lower_bound(inequality_constraint_index));
            // damping
            if (!upper_bounded) {
               barrier_terms += this->damping_factor*(iterate.primals[slack_index] - this->slack_lower_bound(inequality_constraint_index));
            }
         }
         if (upper_bounded) {
            barrier_terms -= std::log(this->slack_upper_bound(inequality_constraint_index) - iterate.primals[slack_index]);
            // damping
            if (!lower_bounded) {
               barrier_terms += this->damping_factor*(this->slack_upper_bound(inequality_constraint_index) - iterate.primals[slack_index]);
            }
         }
         slack_index++;
      }
      barrier_terms *= this->barrier_parameter;
      assert(!std::isnan(barrier_terms) && "The auxiliary measure is not an number.");
      return barrier_terms;
   }

   double PrimalDualInteriorPointProblem::variable_lower_bound(size_t /*variable_index*/) const {
      return -INF<double>;
   }

   double PrimalDualInteriorPointProblem::variable_upper_bound(size_t /*variable_index*/) const {
      return INF<double>;
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_lower_bounded_variables() const {
      return this->first_reformulation.get_lower_bounded_variables();
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_upper_bounded_variables() const {
      return this->first_reformulation.get_upper_bounded_variables();
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_single_lower_bounded_variables() const {
      return this->first_reformulation.get_single_lower_bounded_variables();
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_single_upper_bounded_variables() const {
      return this->first_reformulation.get_single_upper_bounded_variables();
   }

   const Vector<size_t>& PrimalDualInteriorPointProblem::get_fixed_variables() const {
      // fixed variables are handled as linear constraints
      return this->no_fixed_variables; // empty set
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

   size_t PrimalDualInteriorPointProblem::number_objective_gradient_nonzeros() const {
      size_t number_nonzeros = this->first_reformulation.number_objective_gradient_nonzeros();
      // barrier contribution of original variables
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         if (is_finite(this->original_variable_lower_bound(variable_index)) || is_finite(this->original_variable_upper_bound(variable_index))) {
            number_nonzeros++;
         }
      }
      // barrier contribution of slack variables
      number_nonzeros += this->number_slack_variables;
      return number_nonzeros;
   }

   size_t PrimalDualInteriorPointProblem::number_jacobian_nonzeros() const {
      // contribution of original variables
      size_t number_nonzeros = this->first_reformulation.number_jacobian_nonzeros();
      // contribution of slack variables
      number_nonzeros += this->number_slack_variables;
      return number_nonzeros;
   }

   size_t PrimalDualInteriorPointProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      size_t number_nonzeros = this->first_reformulation.number_hessian_nonzeros(hessian_model);
      // barrier contribution of original variables
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         if (is_finite(this->original_variable_lower_bound(variable_index)) || is_finite(this->original_variable_upper_bound(variable_index))) {
            number_nonzeros++;
         }
      }
      // barrier contribution of slack variables
      number_nonzeros += this->number_slack_variables;
      return number_nonzeros;
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate,
         const Multipliers& multipliers) const {
      this->first_reformulation.evaluate_lagrangian_gradient(lagrangian_gradient, iterate, multipliers);

      // barrier terms of the original variables
      for (size_t variable_index: Range(this->first_reformulation.number_variables)) {
         const bool lower_bounded = is_finite(this->original_variable_lower_bound(variable_index));
         const bool upper_bounded = is_finite(this->original_variable_upper_bound(variable_index));
         if (lower_bounded || upper_bounded) {
            double barrier_term = 0.;
            if (lower_bounded) {
               barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - this->original_variable_lower_bound(variable_index));
               // damping
               if (!upper_bounded) {
                  barrier_term += this->damping_factor * this->barrier_parameter;
               }
            }
            if (upper_bounded) {
               barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - this->original_variable_upper_bound(variable_index));
               // damping
               if (!lower_bounded) {
                  barrier_term -= this->damping_factor * this->barrier_parameter;
               }
            }
            // the objective contribution of the Lagrangian gradient may be scaled. Barrier terms go into the constraint contribution
            lagrangian_gradient.constraints_contribution[variable_index] += barrier_term;
         }
      }

      // barrier terms of the original variables
      size_t slack_index = this->first_reformulation.number_variables;
      for (size_t inequality_constraint_index: this->first_reformulation.get_inequality_constraints()) {
         double barrier_term = 0.;
         const bool lower_bounded = is_finite(this->slack_lower_bound(inequality_constraint_index));
         const bool upper_bounded = is_finite(this->slack_upper_bound(inequality_constraint_index));
         if (lower_bounded) {
            barrier_term += -this->barrier_parameter/(iterate.primals[slack_index] - this->slack_lower_bound(inequality_constraint_index));
            // damping
            if (!upper_bounded) {
               barrier_term += this->damping_factor * this->barrier_parameter;
            }
         }
         if (upper_bounded) {
            barrier_term += -this->barrier_parameter/(iterate.primals[slack_index] - this->slack_upper_bound(inequality_constraint_index));
            // damping
            if (!lower_bounded) {
               barrier_term -= this->damping_factor * this->barrier_parameter;
            }
         }
         // the objective contribution of the Lagrangian gradient may be scaled. Barrier terms go into the constraint contribution
         lagrangian_gradient.constraints_contribution[slack_index] += barrier_term;
         slack_index++;
      }
   }

   // form the least-square stationarity \nabla f(x) - z
   void PrimalDualInteriorPointProblem::compute_least_square_stationarity(Iterate& iterate, const Multipliers& multipliers,
         Vector<double>& stationarity_residuals) const {
      stationarity_residuals.fill(0.);

      // objective gradient
      SparseVector<double> objective_gradient(this->first_reformulation.number_variables);
      this->first_reformulation.evaluate_objective_gradient(iterate, objective_gradient);
      for (const auto& [variable_index, derivative]: objective_gradient) {
         stationarity_residuals[variable_index] += derivative;
      }

      // bound duals
      for (size_t variable_index: Range(this->number_variables)) {
         stationarity_residuals[variable_index] -= (multipliers.lower_bounds[variable_index] + multipliers.upper_bounds[variable_index]);
      }
   }

   double PrimalDualInteriorPointProblem::complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const {
      return this->first_reformulation.complementarity_error(primals, constraints, multipliers, shift_value, residual_norm);
   }

   double PrimalDualInteriorPointProblem::original_variable_lower_bound(size_t variable_index) const {
      return this->relaxed_lower_bound(this->first_reformulation.variable_lower_bound(variable_index));
   }

   double PrimalDualInteriorPointProblem::original_variable_upper_bound(size_t variable_index) const {
      return this->relaxed_upper_bound(this->first_reformulation.variable_upper_bound(variable_index));
   }

   double PrimalDualInteriorPointProblem::slack_lower_bound(size_t inequality_constraint_index) const {
      return this->relaxed_lower_bound(this->first_reformulation.constraint_lower_bound(inequality_constraint_index));
   }

   double PrimalDualInteriorPointProblem::slack_upper_bound(size_t inequality_constraint_index) const {
      return this->relaxed_upper_bound(this->first_reformulation.constraint_upper_bound(inequality_constraint_index));
   }

   double PrimalDualInteriorPointProblem::relaxed_lower_bound(double lower_bound) const {
      return lower_bound - this->relaxation_factor * std::max(1., std::abs(lower_bound));
   }

   double PrimalDualInteriorPointProblem::relaxed_upper_bound(double upper_bound) const {
      return upper_bound + this->relaxation_factor * std::max(1., std::abs(upper_bound));
   }
} // namespace