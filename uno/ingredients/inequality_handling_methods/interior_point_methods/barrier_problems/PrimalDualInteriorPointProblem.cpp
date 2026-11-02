// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include "PrimalDualInteriorPointProblem.hpp"
#include "ingredients/subproblem/SymbolicNewtonSystem.hpp"
#include "../InteriorPointParameters.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/View.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   static size_t count_fixed_variables(const OptimizationProblem& problem) {
      size_t number_fixed_variables = 0;
      const auto& lower_bounds = problem.get_variables_lower_bounds();
      const auto& upper_bounds = problem.get_variables_upper_bounds();
      for (size_t variable_index: Range(problem.number_variables)) {
         if (lower_bounds[variable_index] == upper_bounds[variable_index]) {
            ++number_fixed_variables;
         }
      }
      return number_fixed_variables;
   }

   PrimalDualInteriorPointProblem::PrimalDualInteriorPointProblem(const OptimizationProblem& problem,
      const InteriorPointParameters& parameters, const double& barrier_parameter, double residual_scaling_threshold):
         OptimizationProblem(problem.model, problem.number_variables + problem.get_inequality_constraints().size(),
         // move the fixed variables to the set of general constraints
            problem.number_constraints + count_fixed_variables(problem), residual_scaling_threshold),
         unslacked_problem(problem),
         barrier_parameter(barrier_parameter),
         parameters(parameters),
         equality_constraints(this->number_constraints),
         barrier_variables_lower_bounds(this->number_variables, -Inf),
         barrier_variables_upper_bounds(this->number_variables, Inf),
         variables_lower_bounds(this->number_variables, -Inf),
         variables_upper_bounds(this->number_variables, Inf),
         constraints_lower_bounds(this->number_constraints, 0.),
         constraints_upper_bounds(this->number_constraints, 0.),
         constraints_buffer(this->number_constraints),
         predicted_reduction_constraints(this->number_constraints),
         predicted_reduction_jacobian_direction(this->number_constraints) {
      // copy the variables bounds
      const size_t dimension = this->unslacked_problem.number_variables;
      view(this->variables_lower_bounds, 0, dimension) = this->unslacked_problem.get_variables_lower_bounds();
      view(this->variables_upper_bounds, 0, dimension) = this->unslacked_problem.get_variables_upper_bounds();

      // register the inequality constraint of each slack
      size_t inequality_index = 0;
      for (const size_t constraint_index: this->unslacked_problem.get_inequality_constraints()) {
         const size_t slack_index = this->unslacked_problem.number_variables + inequality_index;
         this->slacks.insert(constraint_index, slack_index);
         this->variables_lower_bounds[slack_index] = this->unslacked_problem.get_constraints_lower_bounds()[constraint_index];
         this->variables_upper_bounds[slack_index] = this->unslacked_problem.get_constraints_upper_bounds()[constraint_index];
         ++inequality_index;
      }
      // construct the list of equality constraints (first the original equality constraints)
      std::iota(this->equality_constraints.begin(), this->equality_constraints.begin() +
         static_cast<std::vector<size_t>::difference_type>(problem.number_constraints), /* start */ 0);

      // handle the fixed variables
      size_t fixed_variable_constraint_index = problem.number_constraints;
      size_t equality_constraint_index = problem.number_constraints;
      const auto& lower_bounds = problem.get_variables_lower_bounds();
      const auto& upper_bounds = problem.get_variables_upper_bounds();
      this->fixed_variables.reserve(count_fixed_variables(problem));
      for (size_t variable_index: Range(problem.number_variables)) {
         if (lower_bounds[variable_index] == upper_bounds[variable_index]) { // fixed variable
            // relax the bounds of the fixed variables
            this->variables_lower_bounds[variable_index] = -Inf;
            this->variables_upper_bounds[variable_index] = Inf;
            this->equality_constraints[equality_constraint_index] = fixed_variable_constraint_index;
            this->fixed_variables.insert(fixed_variable_constraint_index, variable_index);
            ++fixed_variable_constraint_index;
            ++equality_constraint_index;
         }
      }

      // variables with barrier terms
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index]) || is_finite(this->variables_upper_bounds[variable_index])) {
            this->bounded_variables.push_back(variable_index);
         }
      }
      for (const size_t variable_index: this->bounded_variables) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            this->lower_bounded_variables.push_back(variable_index);
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            this->upper_bounded_variables.push_back(variable_index);
         }
      }
      this->bound_dual_direction_buffer.resize(this->lower_bounded_variables.size() + this->upper_bounded_variables.size());

      // compute the Jacobian sparsity
      const size_t number_jacobian_nonzeros = this->unslacked_problem.number_jacobian_nonzeros();
      this->jacobian_row_indices.resize(number_jacobian_nonzeros + this->slacks.size() + this->fixed_variables.size());
      this->jacobian_column_indices.resize(number_jacobian_nonzeros + this->slacks.size() + this->fixed_variables.size());
      // first copy the Jacobian of the unslacked problem
      view(this->jacobian_row_indices, 0, number_jacobian_nonzeros) = this->unslacked_problem.get_jacobian_row_indices();
      view(this->jacobian_column_indices, 0, number_jacobian_nonzeros) = this->unslacked_problem.get_jacobian_column_indices();
      size_t nonzero_index = number_jacobian_nonzeros;
      // then the slacks
      for (const auto [constraint_index, slack_index]: this->slacks) {
         this->jacobian_row_indices[nonzero_index] = static_cast<uno_int>(constraint_index);
         this->jacobian_column_indices[nonzero_index] = static_cast<uno_int>(slack_index);
         ++nonzero_index;
      }
      // then the fixed variables
      for (const auto [constraint_index, variable_index]: this->fixed_variables) {
         this->jacobian_row_indices[nonzero_index] = static_cast<uno_int>(constraint_index); // constraint
         this->jacobian_column_indices[nonzero_index] = static_cast<uno_int>(variable_index); // variable
         ++nonzero_index;
      }
   }

   double PrimalDualInteriorPointProblem::get_objective_multiplier() const {
      return this->unslacked_problem.get_objective_multiplier();
   }

   bool PrimalDualInteriorPointProblem::has_inequality_constraints() const {
      return false;
   }

   bool PrimalDualInteriorPointProblem::has_bound_constraints() const {
      return false;
   }

   void PrimalDualInteriorPointProblem::create_initial_iterate(Iterate& iterate, Evaluations& evaluations) const {
      iterate.set_dimensions(this->number_variables, this->number_constraints);

      // make the initial point strictly feasible wrt the bounds
      bool iterate_changed = false;
      for (size_t variable_index: Range(this->unslacked_problem.number_variables)) {
         const double old_value = iterate.primals[variable_index];
         iterate.primals[variable_index] = push_variable_to_interior(iterate.primals[variable_index],
            this->variables_lower_bounds[variable_index], this->variables_upper_bounds[variable_index],
            this->parameters.push_variable_to_interior_k1, this->parameters.push_variable_to_interior_k2);
         if (iterate.primals[variable_index] != old_value) {
            iterate_changed = true;
         }
      }
      // fixed variables
      const auto& fixed_values = this->unslacked_problem.get_variables_lower_bounds();
      for (const auto [_, variable_index]: this->fixed_variables) {
         const double old_value = iterate.primals[variable_index];
         iterate.primals[variable_index] = fixed_values[variable_index];
         if (iterate.primals[variable_index] != old_value) {
            iterate_changed = true;
         }
      }
      if (iterate_changed) {
         evaluations.reset();
      }

      // set the slack variables (if any)
      if (!this->slacks.is_empty()) {
         evaluations.evaluate_constraints(this->model, iterate.primals);
         // set the slacks to the constraint values
         for (const auto [constraint_index, slack_index]: this->slacks) {
            iterate.primals[slack_index] = push_variable_to_interior(evaluations.constraints[constraint_index],
               this->variables_lower_bounds[slack_index], this->variables_upper_bounds[slack_index],
               this->parameters.push_variable_to_interior_k1, this->parameters.push_variable_to_interior_k2);
         }
      }

      // set the bound multipliers
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            iterate.multipliers.lower_bounds[variable_index] = this->parameters.default_multiplier;
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            iterate.multipliers.upper_bounds[variable_index] = -this->parameters.default_multiplier;
         }
      }
   }

   void PrimalDualInteriorPointProblem::postprocess_solution(Iterate& iterate, Evaluations& evaluations) const {
      // fixed variables
      const auto& fixed_bounds = this->unslacked_problem.get_variables_lower_bounds();
      bool iterate_changed = false;
      for (const auto [constraint_index, variable_index]: this->fixed_variables) {
         // fix the primal variable (protects against roundoff)
         if (iterate.primals[variable_index] != fixed_bounds[variable_index]) {
            iterate.primals[variable_index] = fixed_bounds[variable_index];
            iterate_changed = true;
         }
         // set the bound dual
         const double multiplier = iterate.multipliers.constraints[constraint_index];
         iterate.multipliers.lower_bounds[variable_index] = std::max(multiplier, 0.);
         iterate.multipliers.upper_bounds[variable_index] = std::min(multiplier, 0.);
      }
      if (iterate_changed) {
         evaluations.reset();
      }
      // remove the slacks and additional constraints
      iterate.set_dimensions(this->unslacked_problem.number_variables, this->unslacked_problem.number_constraints);

      this->unslacked_problem.postprocess_solution(iterate, evaluations);
   }

   size_t PrimalDualInteriorPointProblem::number_jacobian_nonzeros() const {
      return this->unslacked_problem.number_jacobian_nonzeros() + this->slacks.size() + this->fixed_variables.size();
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
      size_t number_nonzeros = this->unslacked_problem.number_hessian_nonzeros(hessian_model);
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

   void PrimalDualInteriorPointProblem::compute_hessian_sparsity(const HessianModel& hessian_model, View<uno_int> row_indices,
         View<uno_int> column_indices, uno_int solver_indexing) const {
      // original Lagrangian Hessian
      this->unslacked_problem.compute_hessian_sparsity(hessian_model, row_indices, column_indices, solver_indexing);

      // diagonal barrier terms
      size_t current_index = this->unslacked_problem.number_hessian_nonzeros(hessian_model);
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index]) || is_finite(this->variables_upper_bounds[variable_index])) {
            row_indices[current_index] = static_cast<uno_int>(variable_index) + solver_indexing;
            column_indices[current_index] = static_cast<uno_int>(variable_index) + solver_indexing;
            ++current_index;
         }
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_constraints(const Iterate& iterate, View<double> constraints, Evaluations& evaluations) const {
      this->unslacked_problem.evaluate_constraints(iterate, constraints, evaluations);

      // inequality constraints: add the slacks. This makes them homogeneous (c(x) - s = 0)
      for (const auto [constraint_index, slack_index]: this->slacks) {
         constraints[constraint_index] -= iterate.primals[slack_index];
      }

      // make sure the equality constraints are homogeneous (c(x) = 0)
      for (const size_t constraint_index: this->unslacked_problem.get_equality_constraints()) {
         const double fixed_bound = this->unslacked_problem.get_constraints_lower_bounds()[constraint_index];
         constraints[constraint_index] -= fixed_bound;
      }

      // add the fixed variables
      const auto& fixed_bounds = this->unslacked_problem.get_variables_lower_bounds();
      for (const auto [constraint_index, variable_index]: this->fixed_variables) {
         constraints[constraint_index] = iterate.primals[variable_index] - fixed_bounds[variable_index];
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_objective_gradient(const Iterate& iterate, View<double> objective_gradient,
         Evaluations& evaluations) const {
      this->unslacked_problem.evaluate_objective_gradient(iterate, objective_gradient, evaluations);

      // barrier terms + damping
      for (size_t variable_index: Range(this->number_variables)) {
         double barrier_term = 0.;
         if (is_finite(this->variables_lower_bounds[variable_index])) { // lower bounded
            barrier_term -= this->barrier_parameter/(iterate.primals[variable_index] - this->variables_lower_bounds[variable_index]);
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) { // upper bounded
            barrier_term -= this->barrier_parameter/(iterate.primals[variable_index] - this->variables_upper_bounds[variable_index]);
         }
         barrier_term += this->damping_term(variable_index);
         objective_gradient[variable_index] += barrier_term;
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_jacobian(const Vector<double>& primals, View<double> jacobian_values,
         Evaluations& evaluations) const {
      this->unslacked_problem.evaluate_jacobian(primals, jacobian_values, evaluations);

      // add the slack contributions
      size_t nonzero_index = this->unslacked_problem.number_jacobian_nonzeros();
      for ([[maybe_unused]] const auto _: this->slacks) {
         jacobian_values[nonzero_index] = -1.;
         ++nonzero_index;
      }

      // add the contributions of the fixed variables
      for ([[maybe_unused]] const auto _: this->fixed_variables) {
         jacobian_values[nonzero_index] = 1.;
         ++nonzero_index;
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_gradient(const Vector<double>& primals, const Vector<double>& constraint_multipliers,
         const Vector<double>& lower_bound_multipliers, const Vector<double>& upper_bound_multipliers, Evaluations& evaluations,
         View<double> lagrangian_gradient) const {
      this->unslacked_problem.evaluate_lagrangian_gradient(primals, constraint_multipliers, lower_bound_multipliers,
         upper_bound_multipliers, evaluations, lagrangian_gradient);

      for (const auto [_, slack_index]: this->slacks) {
         lagrangian_gradient[slack_index] = 0.;
      }

      // slacks
      for (const auto [constraint_index, slack_index]: this->slacks) {
         // bound multipliers
         lagrangian_gradient[slack_index] -= (lower_bound_multipliers[slack_index] + upper_bound_multipliers[slack_index]);
         // Jacobian block
         lagrangian_gradient[slack_index] += constraint_multipliers[constraint_index];
      }

      // fixed variables
      for (const auto [constraint_index, variable_index]: this->fixed_variables) {
         lagrangian_gradient[variable_index] -= constraint_multipliers[constraint_index];
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, View<double> hessian_values) const {
      // original Lagrangian Hessian
      this->unslacked_problem.evaluate_lagrangian_hessian(statistics, hessian_model, primal_variables, multipliers, hessian_values);

      // barrier terms
      size_t nonzero_index = this->unslacked_problem.number_hessian_nonzeros(hessian_model);
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

   void PrimalDualInteriorPointProblem::add_jacobian_vector_product(View<const double> vector, View<double> result,
         const Evaluations& evaluations) const {
      this->unslacked_problem.add_jacobian_vector_product(vector, result, evaluations);

      // add the slack contributions
      for (const auto [constraint_index, slack_index]: this->slacks) {
         result[constraint_index] -= vector[slack_index];
      }

      // add the contributions of the fixed variables
      for (const auto [constraint_index, variable_index]: this->fixed_variables) {
         result[constraint_index] += vector[variable_index];
      }
   }

   void PrimalDualInteriorPointProblem::add_jacobian_transposed_vector_product(View<const double> vector, View<double> result,
         const Evaluations& evaluations) const {
      this->unslacked_problem.add_jacobian_transposed_vector_product(vector, result, evaluations);

      // add the slack contributions
      for (const auto [constraint_index, slack_index]: this->slacks) {
         result[slack_index] -= vector[constraint_index];
      }

      // add the contributions of the fixed variables
      for (const auto [constraint_index, variable_index]: this->fixed_variables) {
         result[variable_index] += vector[constraint_index];
      }
   }

   void PrimalDualInteriorPointProblem::compute_hessian_vector_product(HessianModel& hessian_model, View<const double> x,
         View<const double> vector, const Multipliers& multipliers, View<double> result) const {
      // original Lagrangian Hessian
      this->unslacked_problem.compute_hessian_vector_product(hessian_model, x, vector, multipliers, result);

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
      return this->barrier_variables_lower_bounds;
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_variables_upper_bounds() const {
      return this->barrier_variables_upper_bounds;
   }

   size_t PrimalDualInteriorPointProblem::get_number_bounded_variables() const {
      size_t number_bounded_variables = 0;
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            ++number_bounded_variables;
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            ++number_bounded_variables;
         }
      }
      return number_bounded_variables;
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_constraints_lower_bounds() const {
      return this->constraints_lower_bounds;
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_constraints_upper_bounds() const {
      return this->constraints_upper_bounds;
   }

   const std::vector<size_t>& PrimalDualInteriorPointProblem::get_equality_constraints() const {
      return this->equality_constraints;
   }

   const std::vector<size_t>& PrimalDualInteriorPointProblem::get_inequality_constraints() const {
      return this->inequality_constraints;
   }

   Inertia PrimalDualInteriorPointProblem::get_inertia() const {
      return {this->number_variables, this->number_constraints, 0};
   }

   // linearized complementarity of the lower and upper bounds: Z Δx + S Δz = μe - S z with S = X - X_B, i.e. for each
   // bound type the blocks Z and S in the rows of the bound duals, the constant μe of the right-hand side and the current
   // values z (the right-hand side μe - S z is in residual form)
   auto PrimalDualInteriorPointProblem::complementarity_blocks(const Vector<double>& primals, const Multipliers& multipliers) const {
      // the iterates are strictly interior: the distances to the bounds S_L = X - X_L are positive and S_U = X - X_U negative
      // (sign preconditions of the pivots of the condensation)
      const auto blocks = [&](auto distance_sign, const std::vector<size_t>& bounded_variables, const std::vector<double>& bounds,
            const Vector<double>& bound_duals) {
         const size_t number_bounds = bounded_variables.size();
         const auto variable = [&bounded_variables](size_t index) { return bounded_variables[index]; };
         const auto distance_to_bound = [&bounded_variables, &bounds, &primals](size_t index) {
            const size_t variable_index = bounded_variables[index];
            return primals[variable_index] - bounds[variable_index];
         };
         const auto bound_dual = [&bounded_variables, &bound_duals](size_t index) { return bound_duals[bounded_variables[index]]; };
         return std::make_tuple(
            partial_permutation(number_bounds, this->number_variables, number_bounds, IdentityIndex{}, variable, bound_dual), // Z
            diagonal<decltype(distance_sign)::value>(number_bounds, distance_to_bound), // S
            component_vector(number_bounds, [this](size_t /*index*/) { return this->barrier_parameter; }), // μe
            bound_dual); // z
      };
      return std::make_tuple(
         blocks(std::integral_constant<EntrySign, EntrySign::positive>{}, this->lower_bounded_variables, this->variables_lower_bounds,
            multipliers.lower_bounds),
         blocks(std::integral_constant<EntrySign, EntrySign::negative>{}, this->upper_bounded_variables, this->variables_upper_bounds,
            multipliers.upper_bounds));
   }

   // Newton system of the barrier problem in (Δx, -Δy, Δz_L, Δz_U):
   //    [∇²L + δ_w I   J^T      -I_L    -I_U ] [  Δx ]   [ -(∇f + damping) + J^T y + z_L + z_U ]  (residual form:)
   //    [J             -δ_c I                ] [ -Δy ] = [ -(c - c_L)                          ]
   //    [Z_L                     S_L         ] [ Δz_L]   [ μe - S_L z_L                        ]
   //    [Z_U                             S_U ] [ Δz_U]   [ μe - S_U z_U                        ]
   // with S_L = X - X_L, S_U = X - X_U, ∇²L the Lagrangian Hessian of the unslacked problem, J the Jacobian of the slacked
   // problem and I_L, I_U the selections of the lower- and upper-bounded variables. The right-hand side is in residual form:
   // the current values z_L, z_U of the bound duals are declared once, and the terms z_L, z_U (= -(-I_L) z_L - (-I_U) z_U)
   // and -S_L z_L, -S_U z_U are generated from the blocks of the matrix.
   // The bound duals are condensed: [∇²L + Σ + δ_w I  *; J  -δ_c I] with Σ = S_L^{-1} Z_L + S_U^{-1} Z_U, and the generated
   // terms cancel symbolically: the condensed right-hand side is -(∇f + damping) + J^T y + μ S_L^{-1} e + μ S_U^{-1} e
   std::unique_ptr<NewtonSystem> PrimalDualInteriorPointProblem::create_newton_system(HessianModel& hessian_model,
         bool use_dual_regularization) const {
      return make_newton_system([this, &hessian_model, use_dual_regularization](const Iterate& iterate) {
         const size_t n = this->number_variables, m = this->number_constraints;
         const size_t number_lower_bounds = this->lower_bounded_variables.size();
         const size_t number_upper_bounds = this->upper_bounded_variables.size();
         // define the blocks
         auto [lower_complementarity, upper_complementarity] = this->complementarity_blocks(iterate.primals, iterate.multipliers);
         auto& [Z_L, S_L, mu_L, z_L] = lower_complementarity;
         auto& [Z_U, S_U, mu_U, z_U] = upper_complementarity;
         const auto lower_variable = [this](size_t index) { return this->lower_bounded_variables[index]; };
         const auto upper_variable = [this](size_t index) { return this->upper_bounded_variables[index]; };
         auto H = LagrangianHessian(this->unslacked_problem, hessian_model, n) + PrimalRegularization(n);
         auto minus_I_L = partial_permutation(n, number_lower_bounds, number_lower_bounds, lower_variable, IdentityIndex{}, ConstantValue{-1.});
         auto minus_I_U = partial_permutation(n, number_upper_bounds, number_upper_bounds, upper_variable, IdentityIndex{}, ConstantValue{-1.});
         auto dual_regularization = DualRegularization(use_dual_regularization ? m : 0);
         // define the block Newton system: matrix, RHS (constant part) and current values of the bound duals
         return kkt_system(
            // clang-format off
            block_matrix({n, m, number_lower_bounds, number_upper_bounds},
               row(std::move(H),               transpose(ConstraintJacobian(*this)), std::move(minus_I_L), std::move(minus_I_U)),
               row(ConstraintJacobian(*this),  std::move(dual_regularization),       zero,                 zero),
               row(std::move(Z_L),             zero,                                 std::move(S_L),       zero),
               row(std::move(Z_U),             zero,                                 zero,                 std::move(S_U))),
            // clang-format on
            block_vector(
               negative_lagrangian_gradient(*this, [this](const Iterate& current_iterate, View<double> objective_gradient,
                     Evaluations& evaluations) {
                  this->evaluate_damped_objective_gradient(current_iterate, objective_gradient, evaluations);
               }),
               NegativeConstraintResidual(*this),
               std::move(mu_L),
               std::move(mu_U)),
            current_values(current_value<2>(z_L), current_value<3>(z_U)));
      }, condense<2, 3>);
   }

   void PrimalDualInteriorPointProblem::assemble_primal_dual_direction(const Iterate& current_iterate, const NewtonSystem& newton_system,
         const Vector<double>& solution, Direction& direction) const {
      // assemble the primal direction and the constraint dual solution
      OptimizationProblem::assemble_primal_dual_direction(current_iterate, newton_system, solution, direction);
      // "fraction-to-boundary" rule for primal variables and bound constraints multipliers
      const double tau = std::max(this->parameters.tau_min, 1. - this->barrier_parameter);
      direction.primal_dual_step_length = primal_fraction_to_boundary(current_iterate.primals, direction.primals, tau);

      // recover the condensed bound duals + fraction-to-boundary
      newton_system.back_substitute(solution.view(), this->bound_dual_direction_buffer.view());
      this->set_bound_dual_direction(direction.multipliers);
      direction.bound_dual_step_length = dual_fraction_to_boundary(current_iterate.multipliers, direction.multipliers, tau);

      DEBUG << "Fraction-to-boundary rules:\n";
      DEBUG << "primal step length = " << direction.primal_dual_step_length << '\n';
      DEBUG << "bound dual step length = " << direction.bound_dual_step_length << "\n\n";
   }

   double PrimalDualInteriorPointProblem::dual_regularization_factor() const {
      return std::pow(this->barrier_parameter, this->parameters.dual_regularization_exponent);
   }

   // protected member functions

   double PrimalDualInteriorPointProblem::push_variable_to_interior(double value, double lower_bound, double upper_bound,
         double absolute_factor, double relative_factor) {
      assert(0. < absolute_factor && 0. < relative_factor && relative_factor <= 0.5);
      const bool has_lower_bound = is_finite(lower_bound);
      const bool has_upper_bound = is_finite(upper_bound);
      if (has_lower_bound && has_upper_bound) {
         const double range = upper_bound - lower_bound;
         assert(0. < range && "push_variable_to_interior: the bounds must satisfy lower_bound < upper_bound");
         const double lower_perturbation = std::min(absolute_factor * std::max(1., std::abs(lower_bound)), relative_factor * range);
         const double upper_perturbation = std::min(absolute_factor * std::max(1., std::abs(upper_bound)), relative_factor * range);
         value = std::max(value, lower_bound + lower_perturbation);
         value = std::min(value, upper_bound - upper_perturbation);
      }
      else if (has_lower_bound) {
         value = std::max(value, lower_bound + absolute_factor * std::max(1., std::abs(lower_bound)));
      }
      else if (has_upper_bound) {
         value = std::min(value, upper_bound - absolute_factor * std::max(1., std::abs(upper_bound)));
      }
      return value;
   }

   // IPOPT (IpoptCalculatedQuantities::CalculateSafeSlack): a slack below s_min = eps*min(1, mu) is replaced with
   // min(max(mu/|z|, s_min), max(slack, 0) + slack_move*max(1, |bound|)). The result is always >= s_min > 0
   double PrimalDualInteriorPointProblem::compute_safe_slack(double slack, double bound_multiplier, double barrier_parameter,
         double bound, double slack_move) {
      const double minimum_slack = std::numeric_limits<double>::epsilon() * std::min(1., barrier_parameter);
      if (minimum_slack <= slack) {
         return slack;
      }
      const double centered_slack = (bound_multiplier != 0.) ? barrier_parameter / std::abs(bound_multiplier) : Inf;
      const double maximum_slack = std::max(slack, 0.) + slack_move * std::max(1., std::abs(bound));
      return std::min(std::max(centered_slack, minimum_slack), maximum_slack);
   }

   // IPOPT (OrigIpoptNLP::AdjustVariableBounds): instead of moving the iterate, move the bound so that the slack becomes safe
   void PrimalDualInteriorPointProblem::reset_slack_bounds() const {
      for (const auto [constraint_index, slack_index]: this->slacks) {
         this->variables_lower_bounds[slack_index] = this->unslacked_problem.get_constraints_lower_bounds()[constraint_index];
         this->variables_upper_bounds[slack_index] = this->unslacked_problem.get_constraints_upper_bounds()[constraint_index];
      }
   }

   size_t PrimalDualInteriorPointProblem::relax_bounds_around(const Iterate& iterate) const {
      size_t number_relaxed_bounds = 0;
      for (size_t variable_index: this->bounded_variables) {
         const double primal = iterate.primals[variable_index];
         double& lower_bound = this->variables_lower_bounds[variable_index];
         if (is_finite(lower_bound)) {
            const double slack = primal - lower_bound;
            const double safe_slack = compute_safe_slack(slack, iterate.multipliers.lower_bounds[variable_index],
               this->barrier_parameter, lower_bound, this->parameters.slack_move);
            if (safe_slack != slack) {
               // guard against rounding (|x| >> safe slack): the new bound must lie strictly below x
               const double relaxed_bound = std::min(primal - safe_slack, std::nextafter(primal, -Inf));
               if (relaxed_bound < lower_bound) {
                  DEBUG << "Slack too small: lower bound of variable " << variable_index << " relaxed from " << lower_bound <<
                     " to " << relaxed_bound << '\n';
                  lower_bound = relaxed_bound;
                  ++number_relaxed_bounds;
               }
            }
         }
         double& upper_bound = this->variables_upper_bounds[variable_index];
         if (is_finite(upper_bound)) {
            const double slack = upper_bound - primal;
            const double safe_slack = compute_safe_slack(slack, iterate.multipliers.upper_bounds[variable_index],
               this->barrier_parameter, upper_bound, this->parameters.slack_move);
            if (safe_slack != slack) {
               const double relaxed_bound = std::max(primal + safe_slack, std::nextafter(primal, Inf));
               if (upper_bound < relaxed_bound) {
                  DEBUG << "Slack too small: upper bound of variable " << variable_index << " relaxed from " << upper_bound <<
                     " to " << relaxed_bound << '\n';
                  upper_bound = relaxed_bound;
                  ++number_relaxed_bounds;
               }
            }
         }
      }
      this->number_bound_relaxations += number_relaxed_bounds;
      return number_relaxed_bounds;
   }

   size_t PrimalDualInteriorPointProblem::get_number_bound_relaxations() const {
      return this->number_bound_relaxations;
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_internal_variables_lower_bounds() const {
      return this->variables_lower_bounds;
   }

   const std::vector<double>& PrimalDualInteriorPointProblem::get_internal_variables_upper_bounds() const {
      return this->variables_upper_bounds;
   }

   const SparseVector<size_t>& PrimalDualInteriorPointProblem::get_slacks() const {
      return this->slacks;
   }

   void PrimalDualInteriorPointProblem::postprocess_iterate(Iterate& iterate, Evaluations& /*evaluations*/) const {
      // 1. with tiny slacks (x - lb), x + alpha d may round onto (or past) a bound: relax the bound (the primals and
      // therefore the evaluations are unchanged)
      this->relax_bounds_around(iterate);

      // 2. κ_Σ reset of the bound multipliers (Eq. 16 in Ipopt paper)
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            const double coefficient = this->barrier_parameter / (iterate.primals[variable_index] - this->variables_lower_bounds[variable_index]);
            if (is_finite(coefficient)) {
               const double lb = coefficient / this->parameters.k_sigma;
               const double ub = coefficient * this->parameters.k_sigma;
               assert(lb <= ub && "The bounds of the multiplier reset are in the wrong order (the iterate is not interior)");
               const double current_value = iterate.multipliers.lower_bounds[variable_index];
               iterate.multipliers.lower_bounds[variable_index] = std::max(std::min(iterate.multipliers.lower_bounds[variable_index], ub), lb);
               if (iterate.multipliers.lower_bounds[variable_index] != current_value) {
                  DEBUG3 << "Multiplier for lower bound " << variable_index << " rescaled from " << current_value << " to " <<
                     iterate.multipliers.lower_bounds[variable_index] << '\n';
               }
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            const double coefficient = this->barrier_parameter / (iterate.primals[variable_index] - this->variables_upper_bounds[variable_index]);
            if (is_finite(coefficient)) {
               const double lb = coefficient * this->parameters.k_sigma;
               const double ub = coefficient / this->parameters.k_sigma;
               assert(lb <= ub && "The bounds of the multiplier reset are in the wrong order (the iterate is not interior)");
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

   // max of the constraint violation of the unslacked problem and of the residuals of the slacked constraints c(x) - s
   // (IPOPT measures the latter): a point where c(x) is feasible but the slacks lag behind is not a KKT point, since the
   // complementarity is measured with the slacks
   double PrimalDualInteriorPointProblem::primal_infeasibility(const Iterate& iterate, Evaluations& evaluations,
         Norm residual_norm) const {
      const double unslacked_infeasibility = this->unslacked_problem.primal_infeasibility(iterate, evaluations, residual_norm);
      this->evaluate_constraints(iterate, this->constraints_buffer.view(), evaluations);
      return std::max(unslacked_infeasibility, norm(residual_norm, this->constraints_buffer));
   }

   double PrimalDualInteriorPointProblem::complementarity_error(const Vector<double>& primals, const Vector<double>& /*constraints*/,
         const Multipliers& multipliers, Norm residual_norm) const {
      // dematerialized vectors
      const auto variable_lower_bound_complementarity = [&](size_t variable_index) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            return multipliers.lower_bounds[variable_index] * (primals[variable_index] - this->variables_lower_bounds[variable_index]);
         }
         return 0.;
      };
      const auto variable_upper_bound_complementarity = [&](size_t variable_index) {
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            return multipliers.upper_bounds[variable_index] * (primals[variable_index] - this->variables_upper_bounds[variable_index]);
         }
         return 0.;
      };
      return norm(residual_norm,
         Range(this->number_variables), variable_lower_bound_complementarity,
         Range(this->number_variables), variable_upper_bound_complementarity
      );
   }

   double PrimalDualInteriorPointProblem::compute_stationarity_scaling(const Multipliers& multipliers) const {
      size_t number_lower_bounded_variables = 0;
      size_t number_upper_bounded_variables = 0;
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            ++number_lower_bounded_variables;
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            ++number_upper_bounded_variables;
         }
      }
      const size_t total_size = number_lower_bounded_variables + number_upper_bounded_variables + model.number_constraints;
      if (total_size == 0) {
         return 1.;
      }
      const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
      const double multiplier_norm = norm_1(multipliers.constraints, multipliers.lower_bounds, multipliers.upper_bounds);
      return std::max(1., multiplier_norm / scaling_factor);
   }

   double PrimalDualInteriorPointProblem::compute_complementarity_scaling(const Multipliers& multipliers) const {
      size_t number_lower_bounded_variables = 0;
      size_t number_upper_bounded_variables = 0;
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            ++number_lower_bounded_variables;
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            ++number_upper_bounded_variables;
         }
      }
      const size_t total_size = number_lower_bounded_variables + number_upper_bounded_variables;
      if (total_size == 0) {
         return 1.;
      }
      const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
      const double bound_multiplier_norm = norm_1(multipliers.lower_bounds, multipliers.upper_bounds);
      return std::max(1., bound_multiplier_norm / scaling_factor);
   }

   void PrimalDualInteriorPointProblem::set_infeasibility_measure(Iterate& iterate, Evaluations& evaluations, Norm progress_norm) const {
      this->evaluate_constraints(iterate, this->constraints_buffer.view(), evaluations);
      iterate.progress.infeasibility = norm(progress_norm, this->constraints_buffer);
   }

   void PrimalDualInteriorPointProblem::set_objective_measure(Iterate& iterate, Evaluations& evaluations) const {
      this->unslacked_problem.set_objective_measure(iterate, evaluations);
   }

   void PrimalDualInteriorPointProblem::set_auxiliary_measure(Iterate& iterate) const {
      // start with the auxiliary measure of the initial problem
      this->unslacked_problem.set_auxiliary_measure(iterate);

      if (is_infinite(this->barrier_parameter)) {
         throw std::runtime_error("Barrier parameter is infinite");
      }

      // add the contribution of the barrier terms
      double barrier_terms = 0.;
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            const double distance = iterate.primals[variable_index] - this->variables_lower_bounds[variable_index];
            if (distance <= 0.) {
               iterate.progress.auxiliary = Inf;
               return;
            }
            barrier_terms -= std::log(distance);
            if (is_infinite(this->variables_upper_bounds[variable_index])) {
               // damping
               barrier_terms += this->parameters.damping_factor*distance;
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            const double distance = this->variables_upper_bounds[variable_index] - iterate.primals[variable_index];
            if (distance <= 0.) {
               iterate.progress.auxiliary = Inf;
               return;
            }
            barrier_terms -= std::log(distance);
            if (is_infinite(this->variables_lower_bounds[variable_index])) {
               // damping
               barrier_terms += this->parameters.damping_factor*distance;
            }
         }
      }
      barrier_terms *= this->barrier_parameter;
      if (std::isnan(barrier_terms)) {
         throw std::runtime_error("The auxiliary measure is not an number");
      }
      iterate.progress.auxiliary += barrier_terms;
   }

   // predicted reductions

   PredictedInfeasibilityReduction PrimalDualInteriorPointProblem::build_predicted_infeasibility_reduction(
         const Iterate& current_iterate, const Vector<double>& primal_direction, Norm progress_norm,
         Evaluations& current_evaluations) const {
      // the model reads its own buffers when it is queried: the scratch buffer constraints_buffer is overwritten when the
      // progress measures of the trial iterates are evaluated
      this->evaluate_constraints(current_iterate, this->predicted_reduction_constraints.view(), current_evaluations);
      const double current_constraint_violation = norm(progress_norm, this->predicted_reduction_constraints);

      this->predicted_reduction_jacobian_direction.fill(0.);
      this->add_jacobian_vector_product(primal_direction.view(), this->predicted_reduction_jacobian_direction.view(),
         current_evaluations);

      return [current_constraint_violation, progress_norm, this](double step_length) {
         const auto linearized_constraints = this->predicted_reduction_constraints + step_length *
            this->predicted_reduction_jacobian_direction;
         const double trial_constraint_violation = norm(progress_norm, linearized_constraints);
         return current_constraint_violation - trial_constraint_violation;
      };
   }

   PredictedObjectiveReduction PrimalDualInteriorPointProblem::build_predicted_objective_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction, Evaluations& current_evaluations, double hessian_quadratic_form) const {
      return this->unslacked_problem.build_predicted_objective_reduction(current_iterate, primal_direction,
         current_evaluations, hessian_quadratic_form);
   }

   PredictedAuxiliaryReduction PrimalDualInteriorPointProblem::build_predicted_auxiliary_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction) const {
      // inner auxiliary model + barrier term  α ↦ α(−μ Xˉ¹eᵀd)
      PredictedAuxiliaryReduction inner_model = this->unslacked_problem.build_predicted_auxiliary_reduction(current_iterate, primal_direction);
      const double barrier_directional_derivative = this->compute_barrier_term_directional_derivative(current_iterate, primal_direction);
      return [inner_model = std::move(inner_model), barrier_directional_derivative](double step_length) {
         return inner_model(step_length) + step_length * (-barrier_directional_derivative);
      };
   }

   double PrimalDualInteriorPointProblem::compute_centrality_error(const Vector<double>& primals, const Multipliers& multipliers,
         double shift) const {
      const auto shifted_lower_bound_complementarity = [&](size_t variable_index) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            return multipliers.lower_bounds[variable_index] * (primals[variable_index] - variables_lower_bounds[variable_index]) - shift;
         }
         return 0.;
      };
      const auto shifted_upper_bound_complementarity = [&](size_t variable_index) {
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            return multipliers.upper_bounds[variable_index] * (primals[variable_index] - variables_upper_bounds[variable_index]) - shift;
         }
         return 0.;
      };
      // TODO use a generic norm
      return norm_inf(
         Range(this->number_variables), shifted_lower_bound_complementarity,
         Range(this->number_variables), shifted_upper_bound_complementarity
      );
   }

   // protected member functions

   // linearized complementarity: bound duals for a given primal direction (back substitution of the complementarity rows)
   void PrimalDualInteriorPointProblem::compute_bound_dual_direction(const Vector<double>& current_primals, const Multipliers& current_multipliers,
         const Vector<double>& direction_primals, Multipliers& direction_multipliers, double& bound_dual_step_length) const {
      const auto [lower_complementarity, upper_complementarity] = this->complementarity_blocks(current_primals, current_multipliers);
      const size_t number_lower_bounds = this->lower_bounded_variables.size();
      const auto& [Z_L, S_L, mu_L, z_L] = lower_complementarity;
      const auto& [Z_U, S_U, mu_U, z_U] = upper_complementarity;
      back_substitute_block(Z_L, S_L, mu_L, direction_primals.view(),
         view(this->bound_dual_direction_buffer.data(), number_lower_bounds), z_L);
      back_substitute_block(Z_U, S_U, mu_U, direction_primals.view(),
         view(this->bound_dual_direction_buffer.data() + number_lower_bounds, this->upper_bounded_variables.size()), z_U);
      this->set_bound_dual_direction(direction_multipliers);
      // fraction to boundary
      const double tau = std::max(this->parameters.tau_min, 1. - this->barrier_parameter);
      bound_dual_step_length = dual_fraction_to_boundary(current_multipliers, direction_multipliers, tau);
   }

   // copy the condensed bound duals (lower-bounded, then upper-bounded variables) into the direction and check their finiteness
   void PrimalDualInteriorPointProblem::set_bound_dual_direction(Multipliers& direction_multipliers) const {
      direction_multipliers.lower_bounds.fill(0.);
      direction_multipliers.upper_bounds.fill(0.);
      size_t index = 0;
      for (const size_t variable_index: this->lower_bounded_variables) {
         direction_multipliers.lower_bounds[variable_index] = this->bound_dual_direction_buffer[index++];
         if (is_infinite(direction_multipliers.lower_bounds[variable_index])) {
            throw std::runtime_error("The lower bound dual is infinite");
         }
      }
      for (const size_t variable_index: this->upper_bounded_variables) {
         direction_multipliers.upper_bounds[variable_index] = this->bound_dual_direction_buffer[index++];
         if (is_infinite(direction_multipliers.upper_bounds[variable_index])) {
            throw std::runtime_error("The upper bound dual is infinite");
         }
      }
   }

   // damping of the barrier terms of the variables with a single finite bound
   double PrimalDualInteriorPointProblem::damping_term(size_t variable_index) const {
      const bool finite_lower_bound = is_finite(this->variables_lower_bounds[variable_index]);
      const bool finite_upper_bound = is_finite(this->variables_upper_bounds[variable_index]);
      if (finite_lower_bound && !finite_upper_bound) {
         return this->parameters.damping_factor * this->barrier_parameter;
      }
      else if (!finite_lower_bound && finite_upper_bound) {
         return -this->parameters.damping_factor * this->barrier_parameter;
      }
      return 0.;
   }

   // gradient of the objective of the unslacked problem + damping terms (the logarithmic barrier terms are not included:
   // they appear through the eliminated bound duals in the Newton system)
   void PrimalDualInteriorPointProblem::evaluate_damped_objective_gradient(const Iterate& iterate, View<double> objective_gradient,
         Evaluations& evaluations) const {
      this->unslacked_problem.evaluate_objective_gradient(iterate, objective_gradient, evaluations);
      for (const size_t variable_index: this->bounded_variables) {
         objective_gradient[variable_index] += this->damping_term(variable_index);
      }
   }

   // TODO use a single function for primal and dual fraction-to-boundary rules
   double PrimalDualInteriorPointProblem::primal_fraction_to_boundary(const Vector<double>& current_primals,
         const Vector<double>& primal_direction, double tau) const {
      double step_length = 1.;
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index]) && primal_direction[variable_index] < 0.) {
            const double distance = -tau * (current_primals[variable_index] - this->variables_lower_bounds[variable_index]) /
               primal_direction[variable_index];
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
      assert(0. < step_length && step_length <= 1. && "The fraction-to-boundary step length should be in (0, 1]");
      return step_length;
   }

   double PrimalDualInteriorPointProblem::dual_fraction_to_boundary(const Multipliers& current_multipliers,
         const Multipliers& direction_multipliers, double tau) const {
      double step_length = 1.;
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index]) && direction_multipliers.lower_bounds[variable_index] < 0.) {
            const double distance = -tau * current_multipliers.lower_bounds[variable_index] / direction_multipliers.lower_bounds[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index]) && 0. < direction_multipliers.upper_bounds[variable_index]) {
            const double distance = -tau * current_multipliers.upper_bounds[variable_index] / direction_multipliers.upper_bounds[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      assert(0. < step_length && step_length <= 1. && "The fraction-to-boundary step length should be in (0, 1]");
      return step_length;
   }
   
   double PrimalDualInteriorPointProblem::compute_barrier_term_directional_derivative(const Iterate& current_iterate,
         const Vector<double>& primal_direction) const {
      double directional_derivative = 0.;
      for (size_t variable_index: Range(this->number_variables)) {
         if (is_finite(this->variables_lower_bounds[variable_index])) {
            directional_derivative += -this->barrier_parameter / (current_iterate.primals[variable_index] -
               this->variables_lower_bounds[variable_index]) * primal_direction[variable_index];
            if (is_infinite(this->variables_upper_bounds[variable_index])) {
               // damping
               directional_derivative += this->parameters.damping_factor * this->barrier_parameter * primal_direction[variable_index];
            }
         }
         if (is_finite(this->variables_upper_bounds[variable_index])) {
            directional_derivative += -this->barrier_parameter / (current_iterate.primals[variable_index] -
               this->variables_upper_bounds[variable_index]) * primal_direction[variable_index];
            if (is_infinite(this->variables_lower_bounds[variable_index])) {
               // damping
               directional_derivative -= this->parameters.damping_factor * this->barrier_parameter * primal_direction[variable_index];
            }
         }
      }
      return directional_derivative;
   }
} // namespace
