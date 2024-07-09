// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXEDPROBLEM_H
#define UNO_L1RELAXEDPROBLEM_H

#include "OptimizationProblem.hpp"
#include "ingredients/constraint_relaxation_strategy/ElasticVariables.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/LagrangianGradient.hpp"
#include "symbolic/Expression.hpp"
#include "symbolic/VectorExpression.hpp"
#include "symbolic/Concatenation.hpp"
#include "symbolic/Range.hpp"
#include "tools/Infinity.hpp"

class l1RelaxedProblem: public OptimizationProblem {
public:
   l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient);

   [[nodiscard]] double get_objective_multiplier() const override;
   void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const override;
   void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const override;
   void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const Vector<double>& x, const Vector<double>& multipliers, SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
   [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
   [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override;
   [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override;
   [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override;
   [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override;

   [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
   [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;

   [[nodiscard]] size_t number_objective_gradient_nonzeros() const override;
   [[nodiscard]] size_t number_jacobian_nonzeros() const override;
   [[nodiscard]] size_t number_hessian_nonzeros() const override;

   void evaluate_lagrangian_gradient(Iterate& iterate, const Multipliers& multipliers) const override;
   [[nodiscard]] double complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const override;
   [[nodiscard]] TerminationStatus check_convergence_with_given_tolerance(Iterate& current_iterate, double tolerance) const override;

   // parameterization
   void set_objective_multiplier(double new_objective_multiplier);

   void set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t, double)>& elastic_setting_function) const;

protected:
   double objective_multiplier;
   const double constraint_violation_coefficient;
   ElasticVariables elastic_variables;
   const Concatenation<const Collection<size_t>&, ForwardRange> lower_bounded_variables; // model variables + elastic variables
   const Concatenation<const Collection<size_t>&, ForwardRange> single_lower_bounded_variables; // model variables + elastic variables

   // delegating constructor
   l1RelaxedProblem(const Model& model, ElasticVariables&& elastic_variables, double objective_multiplier, double constraint_violation_coefficient);
};

inline l1RelaxedProblem::l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient):
   // call delegating constructor
   l1RelaxedProblem(model, ElasticVariables::generate(model), objective_multiplier, constraint_violation_coefficient) {
}

// private delegating constructor
inline l1RelaxedProblem::l1RelaxedProblem(const Model& model, ElasticVariables&& elastic_variables, double objective_multiplier,
         double constraint_violation_coefficient):
      OptimizationProblem(model, model.number_variables + elastic_variables.size(), model.number_constraints),
      objective_multiplier(objective_multiplier),
      constraint_violation_coefficient(constraint_violation_coefficient),
      elastic_variables(std::forward<ElasticVariables>(elastic_variables)),
      // lower bounded variables are the model variables + the elastic variables
      lower_bounded_variables(concatenate(this->model.get_lower_bounded_variables(), Range(model.number_variables,
            model.number_variables + this->elastic_variables.size()))),
      single_lower_bounded_variables(concatenate(this->model.get_single_lower_bounded_variables(),
            Range(model.number_variables, model.number_variables + this->elastic_variables.size()))) {
}

inline double l1RelaxedProblem::get_objective_multiplier() const {
   return this->objective_multiplier;
}

inline void l1RelaxedProblem::evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const {
   // scale nabla f(x) by rho
   if (this->objective_multiplier != 0.) {
      iterate.evaluate_objective_gradient(this->model);
      // TODO change this
      objective_gradient = iterate.evaluations.objective_gradient;
      scale(objective_gradient, this->objective_multiplier);
   }
   else {
      objective_gradient.clear();
   }

   // elastic contribution
   for (const auto [_, elastic_index]: this->elastic_variables.positive) {
      objective_gradient.insert(elastic_index, this->constraint_violation_coefficient);
   }
   for (const auto [_, elastic_index]: this->elastic_variables.negative) {
      objective_gradient.insert(elastic_index, this->constraint_violation_coefficient);
   }
}

inline void l1RelaxedProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
   iterate.evaluate_constraints(this->model);
   constraints = iterate.evaluations.constraints;
   // add the contribution of the elastics
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
      constraints[constraint_index] -= iterate.primals[elastic_index];
   }
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
      constraints[constraint_index] += iterate.primals[elastic_index];
   }
}

inline void l1RelaxedProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
   iterate.evaluate_constraint_jacobian(this->model);
   // TODO change this
   constraint_jacobian = iterate.evaluations.constraint_jacobian;
   // add the contribution of the elastics
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
      constraint_jacobian[constraint_index].insert(elastic_index, -1.);
   }
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
      constraint_jacobian[constraint_index].insert(elastic_index, 1.);
   }
}

inline void l1RelaxedProblem::evaluate_lagrangian_hessian(const Vector<double>& x, const Vector<double>& multipliers,
      SymmetricMatrix<double>& hessian) const {
   this->model.evaluate_lagrangian_hessian(x, this->objective_multiplier, multipliers, hessian);

   // extend the dimension of the Hessian by finalizing the remaining columns (note: the elastics do not enter the Hessian)
   for (size_t constraint_index: Range(this->model.number_variables, this->number_variables)) {
      hessian.finalize_column(constraint_index);
   }
}

// Lagrangian gradient split in two parts: objective contribution and constraints' contribution
inline void l1RelaxedProblem::evaluate_lagrangian_gradient(Iterate& iterate, const Multipliers& multipliers) const {
   iterate.lagrangian_gradient.objective_contribution.fill(0.);
   iterate.lagrangian_gradient.constraints_contribution.fill(0.);

   // objective gradient
   for (auto [variable_index, derivative]: iterate.evaluations.objective_gradient) {
      iterate.lagrangian_gradient.objective_contribution[variable_index] += derivative;
   }

   // constraints
   for (size_t constraint_index: Range(iterate.number_constraints)) {
      if (multipliers.constraints[constraint_index] != 0.) {
         for (auto [variable_index, derivative]: iterate.evaluations.constraint_jacobian[constraint_index]) {
            iterate.lagrangian_gradient.constraints_contribution[variable_index] -= multipliers.constraints[constraint_index] * derivative;
         }
      }
   }

   // bound constraints
   for (size_t variable_index: Range(this->model.number_variables)) {
      iterate.lagrangian_gradient.constraints_contribution[variable_index] -= multipliers.lower_bounds[variable_index] + multipliers.upper_bounds[variable_index];
   }
}

// complementary slackness error: expression for violated constraints depends on the definition of the relaxed problem
inline double l1RelaxedProblem::complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
      const Multipliers& multipliers, double shift_value, Norm residual_norm) const {
   // bound constraints
   const VectorExpression bounds_complementarity(Range(this->model.number_variables), [&](size_t variable_index) {
      if (0. < multipliers.lower_bounds[variable_index]) {
         return multipliers.lower_bounds[variable_index] * (primals[variable_index] - this->variable_lower_bound(variable_index)) - shift_value;
      }
      if (multipliers.upper_bounds[variable_index] < 0.) {
         return multipliers.upper_bounds[variable_index] * (primals[variable_index] - this->variable_upper_bound(variable_index)) - shift_value;
      }
      return 0.;
   });

   // general constraints
   const VectorExpression constraints_complementarity(Range(constraints.size()), [&](size_t constraint_index) {
      // lower violated
      if (constraints[constraint_index] < this->constraint_lower_bound(constraint_index)) {
         return (1. - multipliers.constraints[constraint_index]/this->constraint_violation_coefficient) * (constraints[constraint_index] -
            this->constraint_lower_bound(constraint_index)) - shift_value/this->constraint_violation_coefficient;
      }
      // upper violated
      else if (this->constraint_upper_bound(constraint_index) < constraints[constraint_index]) {
         return (1. + multipliers.constraints[constraint_index]/this->constraint_violation_coefficient) * (constraints[constraint_index] -
            this->constraint_upper_bound(constraint_index)) - shift_value/this->constraint_violation_coefficient;
      }
      // inequality
      else if (this->model.get_constraint_bound_type(constraint_index) != EQUAL_BOUNDS) {
         if (0. < multipliers.constraints[constraint_index]) { // lower bound
            return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->constraint_lower_bound(constraint_index)) - shift_value;
         }
         else if (multipliers.constraints[constraint_index] < 0.) { // upper bound
            return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->constraint_upper_bound(constraint_index)) - shift_value;
         }
      }
      return 0.;
   });
   return norm(residual_norm, bounds_complementarity, constraints_complementarity);
}

inline TerminationStatus l1RelaxedProblem::check_convergence_with_given_tolerance(Iterate& current_iterate, double tolerance) const {
   // evaluate termination conditions based on optimality conditions
   const bool stationarity = (current_iterate.residuals.stationarity / current_iterate.residuals.stationarity_scaling <= tolerance);
   const bool complementarity = (current_iterate.residuals.complementarity / current_iterate.residuals.complementarity_scaling <= tolerance);
   const bool primal_feasibility = (current_iterate.residuals.primal_feasibility <= tolerance);

   DEBUG << "\nTermination criteria for l1 relaxed problem with tolerance = " << tolerance << ":\n";
   DEBUG << "Stationarity: " << std::boolalpha << stationarity << '\n';
   DEBUG << "Complementarity: " << std::boolalpha << complementarity << '\n';
   DEBUG << "Primal feasibility: " << std::boolalpha << primal_feasibility << '\n';

   if (this->model.is_constrained() && stationarity && not primal_feasibility && complementarity) {
      // no primal feasibility, stationary point of constraint violation
      return TerminationStatus::INFEASIBLE_STATIONARY_POINT;
   }
   return TerminationStatus::NOT_OPTIMAL;
}

inline double l1RelaxedProblem::variable_lower_bound(size_t variable_index) const {
   if (variable_index < this->model.number_variables) { // model variable
      return this->model.variable_lower_bound(variable_index);
   }
   else { // elastic variable in [0, +inf[
      return 0.;
   }
}

inline double l1RelaxedProblem::variable_upper_bound(size_t variable_index) const {
   if (variable_index < this->model.number_variables) { // model variable
      return this->model.variable_upper_bound(variable_index);
   }
   else { // elastic variable in [0, +inf[
      return INF<double>;
   }
}

inline double l1RelaxedProblem::constraint_lower_bound(size_t constraint_index) const {
   return this->model.constraint_lower_bound(constraint_index);
}

inline double l1RelaxedProblem::constraint_upper_bound(size_t constraint_index) const {
   return this->model.constraint_upper_bound(constraint_index);
}

inline const Collection<size_t>& l1RelaxedProblem::get_lower_bounded_variables() const {
   return this->lower_bounded_variables;
}

inline const Collection<size_t>& l1RelaxedProblem::get_upper_bounded_variables() const {
   // same set as the model
   return this->model.get_upper_bounded_variables();
}

inline const Collection<size_t>& l1RelaxedProblem::get_single_lower_bounded_variables() const {
   return this->single_lower_bounded_variables;
}

inline const Collection<size_t>& l1RelaxedProblem::get_single_upper_bounded_variables() const {
   // same set as the model
   return this->model.get_single_upper_bounded_variables();
}

inline size_t l1RelaxedProblem::number_objective_gradient_nonzeros() const {
   // elastic contribution
   size_t number_nonzeros = this->elastic_variables.size();

   // objective contribution
   if (this->objective_multiplier != 0.) {
      number_nonzeros += this->model.number_objective_gradient_nonzeros();
   }
   return number_nonzeros;
}

inline size_t l1RelaxedProblem::number_jacobian_nonzeros() const {
   return this->model.number_jacobian_nonzeros() + this->elastic_variables.size();
}

inline size_t l1RelaxedProblem::number_hessian_nonzeros() const {
   return this->model.number_hessian_nonzeros();
}

inline void l1RelaxedProblem::set_objective_multiplier(double new_objective_multiplier) {
   assert(0. <= new_objective_multiplier && "The objective multiplier should be non-negative");
   this->objective_multiplier = new_objective_multiplier;
}

inline void l1RelaxedProblem::set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t, double)>&
      elastic_setting_function) const {
   iterate.set_number_variables(this->number_variables);
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
      elastic_setting_function(iterate, constraint_index, elastic_index, -1.);
   }
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
      elastic_setting_function(iterate, constraint_index, elastic_index, 1.);
   }
}

#endif // UNO_L1RELAXEDPROBLEM_H
