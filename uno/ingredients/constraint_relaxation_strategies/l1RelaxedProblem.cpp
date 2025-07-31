// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "l1RelaxedProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/LagrangianGradient.hpp"
#include "symbolic/Concatenation.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   l1RelaxedProblem::l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient,
            double proximal_coefficient, double const* proximal_center):
         OptimizationProblem(model, model.number_variables + model.get_inequality_constraints().size() +
            2*model.get_equality_constraints().size(), model.number_constraints),
         number_elastic_variables(model.get_inequality_constraints().size() + 2*model.get_equality_constraints().size()),
         objective_multiplier(objective_multiplier),
         constraint_violation_coefficient(constraint_violation_coefficient),
         proximal_coefficient(proximal_coefficient),
         proximal_center(proximal_center),
         // lower bounded variables are the model variables + the elastic variables
         lower_bounded_variables(concatenate(this->model.get_lower_bounded_variables(), Range(model.number_variables,
               this->number_variables))),
         single_lower_bounded_variables(concatenate(this->model.get_single_lower_bounded_variables(),
               Range(model.number_variables, this->number_variables))) {
   }

   l1RelaxedProblem::l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient):
      l1RelaxedProblem(model, objective_multiplier, constraint_violation_coefficient, 0., nullptr) {
   }

   double l1RelaxedProblem::get_objective_multiplier() const {
      return this->objective_multiplier;
   }

   void l1RelaxedProblem::evaluate_objective_gradient(Iterate& iterate, Vector<double>& objective_gradient) const {
      // scale nabla f(x) by rho
      if (this->objective_multiplier != 0.) {
         iterate.evaluate_objective_gradient(this->model);
         // TODO change this
         objective_gradient = iterate.evaluations.objective_gradient;
         objective_gradient.scale(this->objective_multiplier);
      }
      else {
         objective_gradient.fill(0.);
      }

      // constraint violation (through elastic variables) contribution
      for (size_t elastic_index: Range(this->model.number_variables, this->number_variables)) {
         objective_gradient[elastic_index] = this->constraint_violation_coefficient;
      }

      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double proximal_term = this->proximal_coefficient * scaling * scaling * (iterate.primals[variable_index] - this->proximal_center[variable_index]);
            objective_gradient[variable_index] += proximal_term;
         }
      }
   }

   void l1RelaxedProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
      iterate.evaluate_constraints(this->model);
      constraints = iterate.evaluations.constraints;

      // add the contribution of the elastic variables
      size_t elastic_index = this->model.number_variables;
      for (size_t inequality_index: this->model.get_inequality_constraints()) {
         if (is_finite(this->model.constraint_lower_bound(inequality_index))) { // negative part
            constraints[inequality_index] += iterate.primals[elastic_index];
         }
         else { // positive part
            constraints[inequality_index] -= iterate.primals[elastic_index];
         }
         elastic_index++;
      }
      for (size_t equality_index: this->model.get_equality_constraints()) {
         constraints[equality_index] += (iterate.primals[elastic_index] - iterate.primals[elastic_index+1]);
         elastic_index += 2;
      }
   }

   void l1RelaxedProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
      iterate.evaluate_constraint_jacobian(this->model);
      // TODO change this
      constraint_jacobian = iterate.evaluations.constraint_jacobian;

      // add the contribution of the elastic variables
      size_t elastic_index = this->model.number_variables;
      for (size_t inequality_index: this->model.get_inequality_constraints()) {
         if (is_finite(this->model.constraint_lower_bound(inequality_index))) { // negative part
            constraint_jacobian[inequality_index].insert(elastic_index, 1.);
         }
         else { // positive part
            constraint_jacobian[inequality_index].insert(elastic_index, -1.);
         }
         elastic_index++;
      }
      for (size_t equality_index: this->model.get_equality_constraints()) {
         constraint_jacobian[equality_index].insert(elastic_index, 1.);
         constraint_jacobian[equality_index].insert(elastic_index+1, -1.);
         elastic_index += 2;
      }
   }

   // Lagrangian gradient split in two parts: objective contribution and constraints' contribution
   void l1RelaxedProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate,
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
      for (size_t variable_index: Range(this->model.number_variables)) {
         lagrangian_gradient.constraints_contribution[variable_index] -= (multipliers.lower_bounds[variable_index] +
            multipliers.upper_bounds[variable_index]);
      }

      // elastic variables
      size_t elastic_index = this->model.number_variables;
      for (size_t constraint_index: Range(this->number_constraints)) {
         if (is_finite(this->model.constraint_lower_bound(constraint_index))) { // negative part
            lagrangian_gradient.constraints_contribution[elastic_index] += this->constraint_violation_coefficient -
               multipliers.constraints[constraint_index] - multipliers.lower_bounds[elastic_index];
            elastic_index++;
         }
         if (is_finite(this->model.constraint_upper_bound(constraint_index))) { // positive part
            lagrangian_gradient.constraints_contribution[elastic_index] += this->constraint_violation_coefficient +
               multipliers.constraints[constraint_index] - multipliers.lower_bounds[elastic_index];
            elastic_index++;
         }
      }

      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double proximal_term = this->proximal_coefficient * scaling * scaling;
            lagrangian_gradient.constraints_contribution[variable_index] += proximal_term * (iterate.primals[variable_index] -
               this->proximal_center[variable_index]);
         }
      }
   }

   void l1RelaxedProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, SymmetricMatrix<size_t, double>& hessian) const {
      hessian_model.evaluate_hessian(statistics, this->model, primal_variables, this->get_objective_multiplier(), multipliers.constraints, hessian);

      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double proximal_term = this->proximal_coefficient * scaling * scaling;
            hessian.insert(variable_index, variable_index, proximal_term);
         }
      }
   }

   void l1RelaxedProblem::compute_hessian_vector_product(HessianModel& hessian_model, const double* vector, const Multipliers& multipliers,
         double* result) const {
      hessian_model.compute_hessian_vector_product(this->model, vector, this->get_objective_multiplier(), multipliers.constraints, result);

      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double proximal_term = this->proximal_coefficient * scaling * scaling;
            result[variable_index] += proximal_term * vector[variable_index];
         }
      }
   }

   IterateStatus l1RelaxedProblem::check_first_order_convergence(const Iterate& current_iterate, double tolerance) const {
      // evaluate termination conditions based on optimality conditions
      const bool feasibility_stationarity = (current_iterate.residuals.stationarity <= tolerance);
      const bool primal_feasibility = (current_iterate.primal_feasibility <= tolerance);
      const bool feasibility_complementarity = (current_iterate.residuals.complementarity <= tolerance);
      const bool no_trivial_duals = current_iterate.feasibility_multipliers.not_all_zero(this->model.number_variables, tolerance);

      DEBUG << "\nTermination criteria for tolerance = " << tolerance << ":\n";
      DEBUG << "Primal feasibility: " << std::boolalpha << primal_feasibility << '\n';
      DEBUG << "Feasibility stationarity: " << std::boolalpha << feasibility_stationarity << '\n';
      DEBUG << "Feasibility complementarity: " << std::boolalpha << feasibility_complementarity << '\n';
      DEBUG << "Not all zero multipliers: " << std::boolalpha << no_trivial_duals << "\n\n";

      if (this->model.is_constrained() && feasibility_stationarity && !primal_feasibility && feasibility_complementarity &&
            no_trivial_duals) {
         // no primal feasibility, stationary point of constraint violation
         return IterateStatus::INFEASIBLE_STATIONARY_POINT;
      }
      return IterateStatus::NOT_OPTIMAL;
   }

   double l1RelaxedProblem::variable_lower_bound(size_t variable_index) const {
      if (variable_index < this->model.number_variables) { // model variable
         return this->model.variable_lower_bound(variable_index);
      }
      else { // elastic variable in [0, +inf[
         return 0.;
      }
   }

   double l1RelaxedProblem::variable_upper_bound(size_t variable_index) const {
      if (variable_index < this->model.number_variables) { // model variable
         return this->model.variable_upper_bound(variable_index);
      }
      else { // elastic variable in [0, +inf[
         return INF<double>;
      }
   }

   const Collection<size_t>& l1RelaxedProblem::get_lower_bounded_variables() const {
      return this->lower_bounded_variables;
   }

   const Collection<size_t>& l1RelaxedProblem::get_upper_bounded_variables() const {
      return this->model.get_upper_bounded_variables();
   }

   const Collection<size_t>& l1RelaxedProblem::get_single_lower_bounded_variables() const {
      return this->single_lower_bounded_variables;
   }

   const Collection<size_t>& l1RelaxedProblem::get_single_upper_bounded_variables() const {
      return this->model.get_single_upper_bounded_variables();
   }

   const Vector<size_t>& l1RelaxedProblem::get_fixed_variables() const {
      return this->model.get_fixed_variables();
   }

   double l1RelaxedProblem::constraint_lower_bound(size_t constraint_index) const {
      return this->model.constraint_lower_bound(constraint_index);
   }

   double l1RelaxedProblem::constraint_upper_bound(size_t constraint_index) const {
      return this->model.constraint_upper_bound(constraint_index);
   }

   const Collection<size_t>& l1RelaxedProblem::get_equality_constraints() const {
      return this->model.get_equality_constraints();
   }

   const Collection<size_t>& l1RelaxedProblem::get_inequality_constraints() const {
      return this->model.get_inequality_constraints();
   }

   const Collection<size_t>& l1RelaxedProblem::get_dual_regularization_constraints() const {
      return this->dual_regularization_constraints;
   }

   size_t l1RelaxedProblem::number_jacobian_nonzeros() const {
      return this->model.number_jacobian_nonzeros() + this->number_elastic_variables;
   }

   size_t l1RelaxedProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      size_t number_nonzeros = hessian_model.number_nonzeros(this->model);
      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         number_nonzeros += this->model.number_variables;
      }
      return number_nonzeros;
   }

   void l1RelaxedProblem::set_proximal_multiplier(double new_proximal_coefficient) {
      this->proximal_coefficient = new_proximal_coefficient;
   }

   void l1RelaxedProblem::set_proximal_center(double const* new_proximal_center) {
      this->proximal_center = new_proximal_center;
   }

   void l1RelaxedProblem::set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t,
         double)>& elastic_setting_function) const {
      iterate.set_number_variables(this->number_variables);
      size_t elastic_index = this->model.number_variables;
      for (size_t inequality_index: this->model.get_inequality_constraints()) {
         if (is_finite(this->model.constraint_lower_bound(inequality_index))) { // negative part
            elastic_setting_function(iterate, inequality_index, elastic_index, 1.);
         }
         else { // positive part
            elastic_setting_function(iterate, inequality_index, elastic_index, -1.);
         }
         elastic_index++;
      }
      for (size_t equality_index: this->model.get_equality_constraints()) {
         elastic_setting_function(iterate, equality_index, elastic_index, 1.);
         elastic_setting_function(iterate, equality_index, elastic_index+1, -1.);
         elastic_index += 2;
      }
   }
} // namespace